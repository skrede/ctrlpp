# Board capture fixtures

Recordings that `scripts/board_capture.sh --self-test` replays, so every refusal
the capture script implements is re-proved on every change without a board.

| File | What it is |
|---|---|
| `evidence.txt` | The boot report of the allocation-evidence image (`ctrlpp_nucleo_h753zi`), captured from a NUCLEO-H753ZI on 2026-10-08 after a mass-storage flash, which decoded as a pin reset. Every check passed. |
| `timing.txt` | The boot report of the timing image (`ctrlpp_nucleo_h753zi_timing`) built from the same sources, captured the same way on the same day. Every check passed. |
| `reference.txt` | The build ids, stack reserve and toolchain of the two images, read out of the ELFs that produced the recordings, and the report lines the recordings predate. The self-test checks against these instead of a built image. |
| `esp32-evidence.txt` | The boot report of the ESP32 leg's allocation-evidence image (`ctrlpp_esp32`), captured from an ESP32-POE-EA on 2026-10-09 under ESP-IDF v6.0.2 after an esptool flash and an EN-pin restart, which the image decoded as `power_on`. Every check passed. |
| `esp32-timing.txt` | The boot report of the ESP32 leg's timing image (`ctrlpp_esp32_timing`), captured the same way in the same session. Every check passed. |
| `esp32-reference.txt` | The ESP32 images' build ids (the sha256 digest of each ELF), the control task's stack as the source names it, and the toolchain, for the two ESP32 recordings. |

Both recordings are the report lines exactly as the board printed them, with
carriage returns removed. The capture header the script writes above a report
was dropped: it carries the station's port, mount and paths, which are not part
of the stream format and which no refusal reads.

There are two recordings per board rather than one because the refusals split
by image: the allocation canary and the stack lines are judged only on the
evidence image, the cycle-counter canary and the timing rows only on the timing
image. The ESP32 leg's one stack line is judged on both.

The self-test never edits these files. It copies the directory to a temporary
root and synthesizes each corruption there, one per refusal status, then
asserts the status and its message. It also asserts that the unchanged
recordings, and three-run series built from them, pass. The ESP32 series get
one synthesized header line, the capture mode a real run's header records,
because the image decodes an EN-pin restart and a power cycle alike as
`power_on`, and `compare` keeps an `en_pin` series apart from a `power_on` one
by that mode alone.

These recordings say nothing about the board on their own. They are a fixed
sample of the stream format. When the firmware's report format changes, a new
green capture replaces the matching file and the references are re-read from
the new images; until then the self-test proves the parser against the format
as it was recorded.

A line the firmware now prints and a recording lacks is named by a `predates=`
key in that recording's reference: `boot_handshake` for the NUCLEO-H753ZI
recordings, which were made before the image waited for the host's byte at
boot, and `heap_self` for the ESP32 recordings, which were made before the heap
view printed its self-count and per-window verdicts. A recording is replayed
without the checks its reference says it predates; a built image predates
nothing. The self-test covers those checks by synthesizing the missing lines
into a copy and checking it as a built image's report, and the key goes when a
new capture replaces the recording.
