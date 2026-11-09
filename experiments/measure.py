#!/usr/bin/env python3
"""Run one command and report its elapsed time and peak memory, host and GPU.

Why this exists instead of ``run-stelar-pro-with-monitor.sh``: that wrapper calls
``/usr/bin/time -v``, which is GNU-only. On macOS ``/usr/bin/time`` is the BSD
build, which rejects ``-v`` outright ("illegal option -- v") and spells the field
"maximum resident set size" only under ``-l`` -- in *bytes*, where GNU reports
*kilobytes*. So the wrapper silently yields no memory figure on darwin and would
mis-scale by 1024 if it ever matched. This script sidesteps both problems, and
``time -v`` could not see GPU memory in any case.

Host memory is measured two independent ways, because neither alone is
trustworthy:

``rusage_peak_bytes``
    ``getrusage(RUSAGE_CHILDREN).ru_maxrss`` after the child is reaped. This is
    the high-water mark of the single largest child process, not a sum. The unit
    is platform-dependent (bytes on darwin, kibibytes on Linux) and is normalized
    here. Exact, but blind to concurrent siblings.

``sampled_peak_bytes``
    The maximum, over samples taken every ``--sample-interval`` seconds, of the
    summed RSS of the whole process subtree. Catches a job whose cost is spread
    over several processes (``run.sh`` plus its JVM), but can miss a spike that
    falls between two samples.

``peak_rss_bytes`` is the larger of the two, which is the host figure to quote.

GPU memory is a separate measurement, not a component of the host one: device
allocations do not appear in RSS at all, so a CUDA run that fills 20 GiB of VRAM
looks the same in ``peak_rss_bytes`` as one that never touches the card. It is
polled in the same loop via ``nvidia-smi``:

``gpu_peak_bytes``
    Peak of the GPU memory attributed to *this job's own process subtree*, from
    ``--query-compute-apps``. Zero while ``gpu_visible`` is true means a GPU was
    there and this job put nothing on it -- which is exactly how a silent
    fallback from CUDA to the CPU backend shows up.

``gpu_device_peak_bytes``
    Peak of the total in use across every visible device, all users included.
    Quote ``gpu_peak_bytes`` as the job's cost, and compare it against this one
    to tell whether the card was shared while measuring.

Both GPU figures require sampling, so ``--sample-interval 0`` disables them
along with ``sampled_peak_bytes``.

Each invocation of this script is a fresh process, so ``RUSAGE_CHILDREN`` covers
exactly one job. Calling ``getrusage`` repeatedly inside a long-lived driver
would instead return a monotonic high-water mark shared by every job that driver
had ever run -- the reason the measurement is isolated in a subprocess at all.

Usage:
    measure.py --json-out run.json [options] -- command [args ...]
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import time

# ru_maxrss is in bytes on darwin and in kibibytes on Linux and the BSDs.
_RU_MAXRSS_TO_BYTES = 1 if sys.platform == "darwin" else 1024

_MIB = 1024 * 1024
_NVIDIA_SMI = shutil.which("nvidia-smi")


def _process_table() -> dict[int, tuple[int, int]]:
    """Map pid -> (ppid, rss_bytes) for every visible process.

    ``ps`` reports RSS in kibibytes on both darwin and Linux. Returns an empty
    map rather than raising if ``ps`` is unavailable, so sampling can never be
    the reason a measured run fails.
    """
    try:
        out = subprocess.run(
            ["ps", "-axo", "pid=,ppid=,rss="],
            capture_output=True, text=True, timeout=15,
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return {}
    table: dict[int, tuple[int, int]] = {}
    for line in out.splitlines():
        fields = line.split()
        if len(fields) < 3:
            continue
        try:
            table[int(fields[0])] = (int(fields[1]), int(fields[2]) * 1024)
        except ValueError:
            continue
    return table


def _subtree_pids(root: int, table: dict[int, tuple[int, int]]) -> set[int]:
    """Every pid in ``root``'s process subtree, ``root`` itself included."""
    children: dict[int, list[int]] = {}
    for pid, (ppid, _) in table.items():
        children.setdefault(ppid, []).append(pid)
    seen: set[int] = set()
    stack = [root]
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        stack.extend(children.get(pid, ()))
    return seen


def _subtree_rss(root: int, table: dict[int, tuple[int, int]]) -> int:
    """Sum the RSS of ``root`` and every descendant present in ``table``."""
    return sum(table[pid][1] for pid in _subtree_pids(root, table) if pid in table)


def _nvidia_smi(fields: str, kind: str) -> list[str] | None:
    """Run one ``nvidia-smi`` CSV query, or ``None`` if the GPU cannot be read.

    ``None`` is reserved for "cannot see the GPU" -- no driver, no binary, a
    failed call. An empty list is a successful query with nothing to report.
    """
    if _NVIDIA_SMI is None:
        return None
    try:
        completed = subprocess.run(
            [_NVIDIA_SMI, f"--query-{kind}={fields}",
             "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=15,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if completed.returncode != 0:
        return None
    return [line for line in completed.stdout.splitlines() if line.strip()]


def _gpu_compute_apps() -> dict[int, int] | None:
    """Map pid -> GPU memory held in bytes, for every CUDA context on the box.

    Three outcomes are deliberately distinguishable: ``None`` (no visible GPU),
    an empty map (a GPU with nothing running on it), and a map holding no pid of
    ours (someone else's job, or a container whose pid namespace differs from
    the one ``ps`` reports).
    """
    lines = _nvidia_smi("pid,used_gpu_memory", "compute-apps")
    if lines is None:
        return None
    apps: dict[int, int] = {}
    for line in lines:
        fields = [field.strip() for field in line.split(",")]
        if len(fields) < 2:
            continue
        try:
            apps[int(fields[0])] = int(fields[1]) * _MIB
        except ValueError:
            continue    # "[N/A]", as reported under MPS and in some containers
    return apps


def _gpu_device_used() -> tuple[int, int] | None:
    """Memory in use across all visible devices: ``(bytes, device_count)``."""
    lines = _nvidia_smi("memory.used", "gpu")
    if lines is None:
        return None
    total = 0
    count = 0
    for line in lines:
        try:
            total += int(line.strip()) * _MIB
        except ValueError:
            continue
        count += 1
    return total, count


class _Sampler:
    """Poll the subtree's host RSS and its GPU memory until told to stop.

    This deliberately *owns* a thread rather than subclassing ``Thread``. A
    subclass shares its namespace with ``Thread``'s own private attributes, and
    ``Thread`` already defines a ``_stop()`` method that ``join()`` calls
    internally -- so an instance attribute named ``_stop`` shadows it and makes
    ``join()`` raise ``TypeError: 'Event' object is not callable``. It surfaces
    only when ``join()`` wins the race to the thread's state lock, so it can
    pass on one machine and fail on another. Composition cannot collide.
    """

    def __init__(self, pid: int, interval: float, trace: bool = False):
        self.pid = pid
        self.interval = interval
        self.peak = 0
        self.samples = 0
        self.gpu_visible = False
        self.gpu_attributed = False
        self.gpu_peak = 0
        self.gpu_device_peak = 0
        self.gpu_device_count = 0
        self.trace: list[tuple[float, int, int, int]] | None = [] if trace else None
        # Give up on the GPU only after repeated failures, and only while it has
        # never answered: one busy nvidia-smi should not blind the whole run.
        self._gpu_blind = _NVIDIA_SMI is None
        self._gpu_failures = 0
        self._origin = time.perf_counter()
        self._halt = threading.Event()
        self._thread = threading.Thread(target=self._loop, daemon=True,
                                        name="rss-sampler")

    def _sample_gpu(self, table: dict[int, tuple[int, int]]) -> tuple[int, int]:
        """One GPU reading: ``(bytes_ours, bytes_on_all_devices)``."""
        apps = _gpu_compute_apps()
        if apps is None:
            self._gpu_failures += 1
            if not self.gpu_visible and self._gpu_failures >= 3:
                self._gpu_blind = True
            return 0, 0
        self.gpu_visible = True
        self._gpu_failures = 0
        mine = _subtree_pids(self.pid, table) & apps.keys()
        if mine:
            self.gpu_attributed = True
        ours = sum(apps[pid] for pid in mine)
        if ours > self.gpu_peak:
            self.gpu_peak = ours
        device = 0
        reading = _gpu_device_used()
        if reading is not None:
            device, self.gpu_device_count = reading
            if device > self.gpu_device_peak:
                self.gpu_device_peak = device
        return ours, device

    def _loop(self) -> None:
        while not self._halt.is_set():
            table = _process_table()
            rss = _subtree_rss(self.pid, table)
            if rss > self.peak:
                self.peak = rss
            ours = device = 0
            if not self._gpu_blind:
                ours, device = self._sample_gpu(table)
            self.samples += 1
            if self.trace is not None:
                self.trace.append(
                    (time.perf_counter() - self._origin, rss, ours, device))
            self._halt.wait(self.interval)

    def start(self) -> None:
        self._thread.start()

    def stop(self) -> None:
        self._halt.set()

    def join(self, timeout: float | None = None) -> None:
        self._thread.join(timeout=timeout)


def _terminate(process: subprocess.Popen, grace: float = 20.0) -> None:
    """Stop the whole job: signal its process group, escalating to SIGKILL.

    The command is launched in its own process group precisely so a timeout can
    reach the JVM, not just the shell that wrapped it.
    """
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(os.getpgid(process.pid), sig)
        except (ProcessLookupError, PermissionError, OSError):
            try:
                process.kill()
            except OSError:
                pass
            return
        try:
            process.wait(timeout=grace)
            return
        except subprocess.TimeoutExpired:
            continue


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run a command, recording elapsed time and peak memory.",
        epilog="Everything after -- is the command to run.",
    )
    parser.add_argument("--json-out", required=True,
                        help="write the measurement record here as JSON")
    parser.add_argument("--stdout", help="redirect the command's stdout to this file")
    parser.add_argument("--stderr", help="redirect the command's stderr to this file")
    parser.add_argument("--cwd", help="working directory for the command")
    parser.add_argument("--label", default="", help="free-form tag copied into the record")
    parser.add_argument("--timeout", type=float,
                        help="kill the job after this many seconds")
    parser.add_argument("--sample-interval", type=float, default=1.0,
                        help="memory sampling period in seconds; 0 disables "
                             "sampling, and with it every GPU figure (default: 1.0)")
    parser.add_argument("--trace",
                        help="also write every sample to this TSV, so the shape "
                             "of memory use over the run can be inspected")
    parser.add_argument("command", nargs=argparse.REMAINDER,
                        help="-- followed by the command and its arguments")

    args = parser.parse_args(argv)
    command = args.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        parser.error("no command given; put it after --")

    out_handle = open(args.stdout, "w", encoding="utf-8") if args.stdout else None
    err_handle = open(args.stderr, "w", encoding="utf-8") if args.stderr else None

    record: dict[str, object] = {
        "label": args.label,
        "command": command,
        "command_display": " ".join(shlex.quote(part) for part in command),
        "cwd": args.cwd or os.getcwd(),
        "platform": sys.platform,
        "nvidia_smi": _NVIDIA_SMI,
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }

    started = time.perf_counter()
    sampler = None
    try:
        process = subprocess.Popen(
            command,
            cwd=args.cwd,
            stdout=out_handle or None,
            stderr=err_handle or None,
            stdin=subprocess.DEVNULL,
            start_new_session=True,      # own process group, so timeouts can reach the JVM
        )
    except OSError as error:
        record.update(exit_code=None, timed_out=False, failed_to_start=str(error),
                      elapsed_seconds=0.0, peak_rss_bytes=0)
        _write(args.json_out, record)
        print(f"measure: could not start command: {error}", file=sys.stderr)
        return 127

    try:
        if args.sample_interval > 0:
            sampler = _Sampler(process.pid, args.sample_interval,
                               trace=bool(args.trace))
            sampler.start()
        timed_out = False
        try:
            exit_code = process.wait(timeout=args.timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            _terminate(process)
            exit_code = process.returncode
    finally:
        # Teardown must never lose the measurement of a run that already
        # finished, so each step is isolated and its failure only noted.
        teardown_error = None
        if sampler is not None:
            try:
                sampler.stop()
                sampler.join(timeout=5.0)
            except BaseException as error:              # noqa: BLE001
                teardown_error = f"{type(error).__name__}: {error}"
        for handle in (out_handle, err_handle):
            if handle is not None:
                try:
                    handle.close()
                except OSError:
                    pass

    elapsed = time.perf_counter() - started
    rusage_peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss \
        * _RU_MAXRSS_TO_BYTES
    sampled_peak = sampler.peak if sampler is not None else 0

    record.update(
        exit_code=exit_code,
        timed_out=timed_out,
        elapsed_seconds=elapsed,
        rusage_peak_bytes=rusage_peak,
        sampled_peak_bytes=sampled_peak,
        peak_rss_bytes=max(rusage_peak, sampled_peak),
        rss_samples=sampler.samples if sampler is not None else 0,
        sample_interval_seconds=args.sample_interval,
        gpu_visible=sampler.gpu_visible if sampler is not None else False,
        gpu_peak_bytes=sampler.gpu_peak if sampler is not None else 0,
        gpu_peak_attributed=sampler.gpu_attributed if sampler is not None else False,
        gpu_device_peak_bytes=sampler.gpu_device_peak if sampler is not None else 0,
        gpu_device_count=sampler.gpu_device_count if sampler is not None else 0,
    )
    if teardown_error is not None:
        record["teardown_error"] = teardown_error
    _write(args.json_out, record)
    if args.trace and sampler is not None and sampler.trace is not None:
        _write_trace(args.trace, sampler.trace)

    summary = (f"measure: elapsed {elapsed:.2f} s, "
               f"peak RSS {record['peak_rss_bytes'] / _MIB:.1f} MiB")
    if sampler is not None and sampler.gpu_visible:
        summary += f", peak GPU {sampler.gpu_peak / _MIB:.1f} MiB"
        if not sampler.gpu_attributed:
            summary += (f" (none attributable to this job; "
                        f"{sampler.gpu_device_peak / _MIB:.1f} MiB device-wide)")
    summary += f", exit {exit_code}{' (TIMED OUT)' if timed_out else ''}"
    print(summary)
    if timed_out:
        return 124
    return exit_code if isinstance(exit_code, int) else 1


def _write(path: str, record: dict[str, object]) -> None:
    directory = os.path.dirname(os.path.abspath(path))
    if directory:
        os.makedirs(directory, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _write_trace(path: str, samples: list[tuple[float, int, int, int]]) -> None:
    """Dump the sample series as TSV. A failure here must not fail the run."""
    try:
        directory = os.path.dirname(os.path.abspath(path))
        if directory:
            os.makedirs(directory, exist_ok=True)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("seconds\trss_mib\tgpu_mib\tgpu_device_mib\n")
            for moment, rss, gpu, device in samples:
                handle.write(f"{moment:.3f}\t{rss / _MIB:.1f}\t"
                             f"{gpu / _MIB:.1f}\t{device / _MIB:.1f}\n")
    except OSError as error:
        print(f"measure: could not write trace {path}: {error}", file=sys.stderr)


if __name__ == "__main__":
    raise SystemExit(main())
