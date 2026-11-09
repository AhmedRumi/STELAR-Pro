#!/usr/bin/env python3
"""Pull per-phase times and memory peaks out of a STELAR-Pro log.

STELAR-Pro already instruments itself more finely in time than any external
wrapper can: ``PhaseLogger`` polls ``/proc/self/status`` VmRSS every 50 ms and,
during GPU phases, ``cudaMemGetInfo()`` every 50 ms, reporting a per-phase and a
run-wide peak. It prints all of that for a human to read and nothing more, so
this module turns it back into rows.

The blocks look like this, one per phase::

      ▶  Phase 1  Parse gene trees  [CPU]
         ✓  8 ms
            peak RAM (this step): 6 MiB  │  peak RAM (so far): 6 MiB

followed at the end of a successful run by ``Total time: N ms`` and a Run
Summary carrying the triplet score and the run-wide peaks.

Two things to keep in mind about the numbers recovered here, both properties of
the source rather than of the parsing:

* ``peak_vram_mib`` is **device-wide**. ``PhaseLogger`` derives it as
  ``total - minFree`` from ``cudaMemGetInfo``, which reports the whole card, so
  another job sharing the GPU inflates it. ``measure.py``'s ``gpu_peak_bytes``
  is the pid-attributed figure; this one has the better time resolution. They
  answer different questions.
* ``peak_ram_mib`` silently changes meaning off Linux. When
  ``/proc/self/status`` cannot be read, ``readVmRssMiB()`` falls back to
  ``totalMemory() - freeMemory()`` -- Java heap only, not RSS, and much smaller.

Usage:
    parse_phases.py run.stelar-pro.log
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

# Colour is normally off when stderr is redirected, but FORCE_COLOR or a tty
# would leave escapes in the file, so they are stripped unconditionally.
_ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")

_BEGIN = re.compile(r"^\s*▶\s+(?P<label>.+?)\s+\[(?P<device>CPU|GPU)\]\s*$")
_END = re.compile(r"^\s*✓\s+(?P<ms>\d+)\s*ms\s*$")
# The separator between the two figures is a box-drawing character; matching
# loosely across it keeps this working if that formatting ever changes.
_PEAK = re.compile(
    r"peak\s+(?P<kind>V?RAM)\s*\(this step\):\s*(?P<step>\d+)\s*MiB"
    r".*?\(so far\):\s*(?P<far>\d+)\s*MiB")
_TOTAL = re.compile(r"Total time:\s*(\d+)\s*ms")
_SUMMARY_SCORE = re.compile(r"^\s*Triplet score\s+(?P<value>\S+)\s*$")
_SUMMARY_RAM = re.compile(r"^\s*Max CPU RAM\s+(?P<value>\d+)\s*MiB\s*$")
_SUMMARY_VRAM = re.compile(r"^\s*Max GPU VRAM\s+(?P<value>\d+)\s*MiB\s*$")


@dataclass
class Phase:
    """One completed phase block."""
    index: int
    label: str
    device: str
    ms: int
    peak_ram_mib: int | None = None
    peak_ram_so_far_mib: int | None = None
    peak_vram_mib: int | None = None
    peak_vram_so_far_mib: int | None = None


@dataclass
class LogSummary:
    """Everything recoverable from one run's log."""
    phases: list[Phase] = field(default_factory=list)
    internal_total_ms: int | None = None
    triplet_score: str | None = None
    log_peak_ram_mib: int | None = None
    log_peak_vram_mib: int | None = None
    unfinished_phase: str | None = None

    @property
    def slowest(self) -> Phase | None:
        return max(self.phases, key=lambda phase: phase.ms, default=None)


def parse_log(path: str | Path) -> LogSummary:
    """Parse a STELAR-Pro log. Never raises: an unreadable log yields an empty
    summary, because a parsing problem must not cost a finished run its row."""
    summary = LogSummary()
    try:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return summary

    open_phase: Phase | None = None
    for raw in text.splitlines():
        line = _ANSI.sub("", raw)

        begun = _BEGIN.match(line)
        if begun:
            # A phase header while one is still open means the previous phase
            # never completed -- worth surfacing rather than quietly dropping.
            if open_phase is not None:
                summary.unfinished_phase = open_phase.label
            open_phase = Phase(index=len(summary.phases) + 1,
                               label=" ".join(begun.group("label").split()),
                               device=begun.group("device"), ms=-1)
            continue

        ended = _END.match(line)
        if ended and open_phase is not None:
            open_phase.ms = int(ended.group("ms"))
            summary.phases.append(open_phase)
            # Cleared here, not left open: the peak lines that follow address
            # phases[-1], so nothing needs this reference, and leaving it set
            # would make the next phase header flag a phase that did finish.
            open_phase = None
            continue

        peak = _PEAK.search(line)
        if peak and summary.phases:
            target = summary.phases[-1]
            if peak.group("kind") == "RAM":
                target.peak_ram_mib = int(peak.group("step"))
                target.peak_ram_so_far_mib = int(peak.group("far"))
            else:
                target.peak_vram_mib = int(peak.group("step"))
                target.peak_vram_so_far_mib = int(peak.group("far"))
            continue

        total = _TOTAL.search(line)
        if total:
            summary.internal_total_ms = int(total.group(1))
            if open_phase is not None and open_phase.ms < 0:
                summary.unfinished_phase = open_phase.label
            open_phase = None
            continue

        for pattern, attribute, cast in (
            (_SUMMARY_SCORE, "triplet_score", str),
            (_SUMMARY_RAM, "log_peak_ram_mib", int),
            (_SUMMARY_VRAM, "log_peak_vram_mib", int),
        ):
            found = pattern.match(line)
            if found:
                setattr(summary, attribute, cast(found.group("value")))
                break

    if open_phase is not None and open_phase.ms < 0:
        summary.unfinished_phase = open_phase.label

    # The Run Summary is only printed for a successful analysis, so fall back to
    # the highest "so far" figure seen in the phase blocks.
    if summary.log_peak_ram_mib is None:
        peaks = [p.peak_ram_so_far_mib for p in summary.phases
                 if p.peak_ram_so_far_mib is not None]
        summary.log_peak_ram_mib = max(peaks) if peaks else None
    if summary.log_peak_vram_mib is None:
        peaks = [p.peak_vram_so_far_mib for p in summary.phases
                 if p.peak_vram_so_far_mib is not None]
        summary.log_peak_vram_mib = max(peaks) if peaks else None
    return summary


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        print(__doc__)
        return 2
    for path in argv:
        summary = parse_log(path)
        print(f"{path}")
        print(f"  internal total   {summary.internal_total_ms} ms")
        print(f"  triplet score    {summary.triplet_score}")
        print(f"  peak RAM         {summary.log_peak_ram_mib} MiB")
        print(f"  peak VRAM        {summary.log_peak_vram_mib} MiB")
        if summary.unfinished_phase:
            print(f"  UNFINISHED       {summary.unfinished_phase}")
        width = max((len(p.label) for p in summary.phases), default=0)
        for phase in summary.phases:
            print(f"    {phase.label:<{width}}  {phase.device}  {phase.ms:>7} ms"
                  f"  RAM {phase.peak_ram_mib} MiB"
                  + (f"  VRAM {phase.peak_vram_mib} MiB"
                     if phase.peak_vram_mib is not None else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
