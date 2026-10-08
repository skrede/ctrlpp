# Board capture fixtures

Recordings that `scripts/board_capture.sh --self-test` replays, so every refusal
the capture script implements is re-proved on every change without a board.

| File | What it is |
|---|---|
| `evidence.txt` | The boot report of the allocation-evidence image (`ctrlpp_nucleo_h753zi`), captured from a NUCLEO-H753ZI on 2026-10-08 after a mass-storage flash, which decoded as a pin reset. Every check passed. |
| `timing.txt` | The boot report of the timing image (`ctrlpp_nucleo_h753zi_timing`) built from the same sources, captured the same way on the same day. Every check passed. |
| `reference.txt` | The build ids, stack reserve and toolchain of the two images, read out of the ELFs that produced the recordings. The self-test checks against these instead of a built image. |

Both recordings are the report lines exactly as the board printed them, with
carriage returns removed. The capture header the script writes above a report
was dropped: it carries the station's port, mount and paths, which are not part
of the stream format and which no refusal reads.

There are two recordings rather than one because the refusals split by image:
the allocation canary and the stack lines are judged only on the evidence image,
the cycle-counter canary and the timing rows only on the timing image.

The self-test never edits these files. It copies the directory to a temporary
root and synthesizes each corruption there, one per refusal status, then
asserts the status and its message. It also asserts that the unchanged
recordings, and three-run series built from them, pass.

These recordings say nothing about the board on their own. They are a fixed
sample of the stream format. When the firmware's report format changes, a new
green capture replaces the matching file and the references are re-read from
the new images; until then the self-test proves the parser against the format
as it was recorded.
