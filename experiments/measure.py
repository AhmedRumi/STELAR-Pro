#!/usr/bin/env python3
"""Run one command and report its elapsed time and peak resident memory.

Why this exists instead of ``run-stelar-pro-with-monitor.sh``: that wrapper calls
``/usr/bin/time -v``, which is GNU-only. On macOS ``/usr/bin/time`` is the BSD
build, which rejects ``-v`` outright ("illegal option -- v") and spells the field
"maximum resident set size" only under ``-l`` -- in *bytes*, where GNU reports
*kilobytes*. So the wrapper silently yields no memory figure on darwin and would
mis-scale by 1024 if it ever matched. This script sidesteps both problems.

Two independent measurements are taken, because neither alone is trustworthy:

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

``peak_rss_bytes`` is the larger of the two, which is the figure to quote.

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
import signal
import subprocess
import sys
import threading
import time

# ru_maxrss is in bytes on darwin and in kibibytes on Linux and the BSDs.
_RU_MAXRSS_TO_BYTES = 1 if sys.platform == "darwin" else 1024


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


def _subtree_rss(root: int, table: dict[int, tuple[int, int]]) -> int:
    """Sum the RSS of ``root`` and every descendant present in ``table``."""
    children: dict[int, list[int]] = {}
    for pid, (ppid, _) in table.items():
        children.setdefault(ppid, []).append(pid)
    total = 0
    stack = [root]
    seen = set()
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        entry = table.get(pid)
        if entry is not None:
            total += entry[1]
        stack.extend(children.get(pid, ()))
    return total


class _Sampler:
    """Poll the subtree's summed RSS until told to stop.

    This deliberately *owns* a thread rather than subclassing ``Thread``. A
    subclass shares its namespace with ``Thread``'s own private attributes, and
    ``Thread`` already defines a ``_stop()`` method that ``join()`` calls
    internally -- so an instance attribute named ``_stop`` shadows it and makes
    ``join()`` raise ``TypeError: 'Event' object is not callable``. It surfaces
    only when ``join()`` wins the race to the thread's state lock, so it can
    pass on one machine and fail on another. Composition cannot collide.
    """

    def __init__(self, pid: int, interval: float):
        self.pid = pid
        self.interval = interval
        self.peak = 0
        self.samples = 0
        self._halt = threading.Event()
        self._thread = threading.Thread(target=self._loop, daemon=True,
                                        name="rss-sampler")

    def _loop(self) -> None:
        while not self._halt.is_set():
            total = _subtree_rss(self.pid, _process_table())
            if total > self.peak:
                self.peak = total
            self.samples += 1
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
        description="Run a command, recording elapsed time and peak RSS.",
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
                        help="RSS sampling period in seconds (default: 1.0)")
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
            sampler = _Sampler(process.pid, args.sample_interval)
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
    )
    if teardown_error is not None:
        record["teardown_error"] = teardown_error
    _write(args.json_out, record)

    peak_mb = record["peak_rss_bytes"] / (1024 * 1024)
    print(f"measure: elapsed {elapsed:.2f} s, peak RSS {peak_mb:.1f} MiB, "
          f"exit {exit_code}{' (TIMED OUT)' if timed_out else ''}")
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


if __name__ == "__main__":
    raise SystemExit(main())
