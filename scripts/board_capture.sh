#!/usr/bin/env bash
#
# board_capture.sh -- build a board's image, then capture one boot report from
# the board and refuse anything that is not evidence.
#
# Usage:
#   board_capture.sh [--board=esp32] prepare
#   board_capture.sh [--board=esp32] capture [--timing] [--reset-only | --no-reset | --power-cycle] <destination-file>
#   board_capture.sh [--board=esp32] check [--timing] <transcript-or-capture-file>
#   board_capture.sh [--board=esp32] compare [--timing --bound=<cycles>] <count> <reset-class> <destination-directory>
#   board_capture.sh [--board=esp32] compare --replay [--timing --bound=<cycles>] <reset-class> <series-directory>
#   board_capture.sh [--board=esp32] discover
#   board_capture.sh --self-test
#
# Boards:
#   (default)      the NUCLEO-H753ZI leg, examples/embedded/nucleo_h753zi, built
#                  here into build-arm-capture. Everything below describes it
#                  unless a paragraph names the ESP32.
#   --board=esp32  the ESP32-POE-EA leg, examples/embedded/esp32, built by the
#                  framework into its build/ directory ("prepare" needs the
#                  framework exported and that tree configured once with
#                  idf.py set-target esp32). Both legs emit the same tagged
#                  two-class stream, so the parse and every refusal are shared;
#                  what differs is listed in "The ESP32 leg" below.
#
# "prepare" configures, builds and produces the raw binaries beside the two
# ELFs, the allocation-evidence image and the timing image. It touches no
# hardware. "check" runs every refusal of the chosen mode over a transcript
# already on disk, against the ELFs the last "prepare" built; it touches no
# hardware either. "capture" performs no build at all: a build discovered
# inside a capture is how a capture ends up describing an image that was never
# flashed, which is the failure this split exists to prevent, so a missing ELF
# or binary is an error naming the target.
#
# "compare" collects a series: one uncounted priming run, then <count> captures of
# the resident image, each restarted by the reset class the series is for, each
# kept as run-NN.txt and each passing every refusal "capture" applies. It then
# compares them and writes compare-summary.txt beside them. --replay compares the
# run-NN.txt files already in a directory, after re-applying "check" to each.
#   reset class    the class the image decodes and prints as reset=<class>. A
#                  'pin' series is restarted with `st-flash reset`. A 'power_on'
#                  series prompts the operator to unplug and replug the board and
#                  waits, bounded, for the serial node to vanish and return; no
#                  tool here can cut the target's power. An 'en_pin' series, on
#                  the ESP32 leg only, is restarted by this script's EN-pin pulse
#                  through the console bridge's RTS line, as capture --reset-only
#                  does, with no prompt; the image decodes that pulse as
#                  reset=power_on, which is what each run must decode, so its
#                  summary and every run's header say the restart was an EN-pin
#                  pulse and not a power cycle, and that the label cannot tell
#                  the two apart. An en_pin series is never power-on evidence and
#                  is never pooled with a power_on series: on the ESP32 leg each
#                  run's header must record the restart its series' class names,
#                  reset-only for en_pin and power-cycle for power_on. A run whose
#                  decoded class is not the required one is not counted: it is
#                  kept as refused-NN.txt, listed in the summary with its decoded
#                  class, and the run is collected again, prompting the operator
#                  where its class needs one, until <count> runs are counted. Any
#                  other refusal of a run ends the series.
#   exact class    the [exact] lines, compared byte for byte in order against the
#                  first run, with no rounding, normalization or sorting; the
#                  refusal names the first differing line and byte.
#   spread class   per [spread] row and field, the minimum, maximum and spread
#                  across the runs. raw_min, raw_median and raw_max are cycle
#                  counts held to --bound; max_at and per_iteration_raw_max are
#                  reported without one; n, n_basis, block and wraps say what was
#                  measured and must agree. A spread is reported even when zero.
#   shared         every run carries one build id and one posture line apart
#                  from its reset fields; the summary quotes both, with the
#                  toolchain, the run count and every run's decoded reset.
#
# "discover" looks for the board the way "capture" does and reports what it
# found, then exits 0 either way: it flashes nothing, reads no report and forms
# no verdict. With no board attached it prints one line beginning "SKIP: no board
# found" that names the on-target part of the check as not performed, so a
# machine without the board says so instead of passing quietly.
#
# "--self-test" replays the committed recordings in board_capture_fixtures/, one
# per image, against the references recorded beside them. It copies the fixture
# tree to a temporary root and, for each refusal status this script implements,
# synthesizes one corruption into the copy and asserts that status and its
# message; it also asserts that the untouched recordings and series built from
# them pass, so the refusals are shown to be selective. It needs no board, no
# port, no flash and no built image.
#
# Images:
#   (default)      the evidence image: carries the allocation sensor, so its
#                  report must carry a passing canary.
#   --timing       the timing image: built from the same sources without the
#                  allocation sensor, so the code it times carries none of the
#                  sentinel's checks. Its report must carry the timing lines and
#                  need not carry a canary; the written capture names the
#                  evidence image's fingerprint beside its own.
#
# Capture modes:
#   (default)      copy the binary onto the ST-LINK mass-storage mount, which
#                  both programs and resets the target.
#   --reset-only   do not program; restart the resident image with plain
#                  `st-flash reset`. This is how a board is re-read without
#                  disturbing what is on it. `--connect-under-reset` is NOT used:
#                  it leaves the core halted and the report never arrives.
#   --no-reset     do not program and do not reset; read whatever the target is
#                  emitting right now. Used when the operator has reset the board
#                  themselves (the B2 button, a power cycle) or when the running
#                  image must not be perturbed.
#   --power-cycle  ESP32 only: do not program; attach the readers, ask the
#                  operator to unplug and replug the board's USB cable, and read
#                  the report of the boot that follows. The telemetry adapter is
#                  a USB device of its own, so its reader is attached through
#                  the cycle and the report keeps its head.
#
# The ESP32 leg:
#   transport      esptool, the framework's own flasher, on the console port
#                  (the board's USB bridge), writes the bootloader, the
#                  partition table and the image named by the build's
#                  flash_args, and verifies what it wrote; the image of the
#                  timing kind takes the application slot in place of the
#                  evidence image. idf.py flash is not used because it builds
#                  first. The board is restarted by its EN pin, driven through
#                  the bridge's RTS line exactly as esptool's own hard reset
#                  does: the readers attach, EN is held low until both ports
#                  are quiet, and only the bytes after that point are read. The
#                  report arrives on a second port, the telemetry adapter, and
#                  the console is recorded beside the capture in
#                  <destination>-console.txt, transported and not judged.
#   identity       the sha256 digest of the ELF, which the framework's build
#                  writes into the image descriptor and the image prints as its
#                  build id; sha256sum of the built ELF reproduces it. Before
#                  flashing, the binary's descriptor is read with esptool and
#                  must name the same digest, so the binary written is the ELF
#                  the report is checked against.
#   report         the evidence image keeps running after its report: the
#                  stack probe that follows it ends in a deliberate overflow and
#                  the panic handler restarts the chip, so that image emits a
#                  report on every boot. The timing image runs no probe and
#                  reports once per boot. The capture is the first report after
#                  the restart this script made, from its build id line to its
#                  closing verdict; residue before that line is discarded and
#                  its size is stated in the capture header. The ESP32 leg has
#                  no boot handshake: its telemetry adapter is a USB device of
#                  its own and stays attached through a restart.
#   differences    binary32 exact rows (8 hexadecimal digits); the predictive
#                  family must be declared ABSENT on exactly one line; no
#                  predictive timing region; posture fields reset_reason and
#                  reset_rom in place of reset_flags and reset_rsr, plus a
#                  schedule line; one stack line, the control task's
#                  high-water, judged against the task's stack.
#
# The reader is attached BEFORE anything resets the target. The image prints its
# report once at boot and then idles, so a reader attached afterwards reads an
# empty stream; that is a failed capture and this script says so rather than
# reporting a pass it did not earn.
#
# Boot handshake (the NUCLEO-H753ZI leg): the image waits at boot, before it
# measures or prints anything, for one byte from the host on the virtual COM
# port, the kBootHandshakeByte that boot_handshake.h names, for at most the
# timeout that header documents. Once its reader is attached, this script sends
# that byte about once a second until the report's build id line arrives, in
# every capture mode, so a reader that attaches late, after a power cycle or a
# re-enumeration, still gets the report's head. The image reports the outcome
# on a [meta] boot_handshake line, and a report whose wait expired is refused:
# it began with no reader attached, so its head may be lost.
#
# The read is bounded twice over: at most CAPTURE_READ_TIMEOUT seconds (default
# 45) and at most CAPTURE_BYTE_CAP bytes (default 262144), whichever comes first,
# and it stops early on the report's closing verdict line.
#
# Stream:
#   Every report line opens with one of five tags, and a block is the set of
#   lines carrying its tag:
#     [meta]     the fingerprint, the image kind, and the posture line, which
#                carries the reset class the image decoded from its own reset
#                status register before clearing it
#     [canary]   the instrument's proof of sight: the allocation sensor on the
#                evidence image, the cycle counter on the timing image
#     [exact]    a header row, then one row per loop index in index order; the
#                computed columns, each double printed as its 16-digit binary64
#                bit pattern, asserted identical across runs
#     [spread]   a header row, then one row per timed region and regime; the
#                timing columns, which carry a spread bound and never a bitwise
#                claim (the evidence image times nothing, so its block is the
#                header alone)
#     [verdict]  every family's verdict, the stack, heap and allocation figures,
#                and the closing golden diff
#   A field belongs to exactly one of the two claim classes, so a name in both
#   headers is refused.
#
# Environment:
#   CAPTURE_BUILD_DIR        build tree (default: <repo>/build-arm-capture)
#   CAPTURE_CONFIGURE_ARGS   NUCLEO-H753ZI only: cache arguments added to the
#                            configure "prepare" runs, for example
#                            -DCTRLPP_MCU_TUNED_POSTURE=ON, which belongs in a
#                            CAPTURE_BUILD_DIR of its own (default: none; a
#                            build tree keeps what it was configured with)
#   CAPTURE_MOUNT            ST-LINK mass-storage mount
#                            (default: /run/media/$USER/NOD_H753ZI)
#   CAPTURE_PORT             serial port, when discovery must be overridden (on
#                            the ESP32 leg, the telemetry port)
#   CAPTURE_CONSOLE_PORT     ESP32 only: the console port, the board's USB
#                            bridge, when discovery must be overridden
#   CAPTURE_CONSOLE_TAIL     ESP32 only: seconds the console is read past the
#                            report's closing line, for the evidence image's
#                            stack probe that prints after it (default: 15)
#   CAPTURE_READ_TIMEOUT     seconds (default: 45)
#   CAPTURE_BYTE_CAP         bytes (default: 262144)
#   CAPTURE_KEEP_RAW         a directory that receives the unparsed streams of
#                            every capture before they are removed, so a refused
#                            capture can be diagnosed (default: unset, discarded)
#   CAPTURE_BUILD_JOBS       build parallelism (default: 6)
#   COMPARE_ATTEMPTS         captures per run when the report head is lost,
#                            status 4 or 5 (default: 3)
#   COMPARE_PLUG_TIMEOUT     seconds the power-on prompt waits for the node to
#                            vanish, and again to return (default: 120)
#   COMPARE_MAX_REFUSED      runs of another class a series may keep before it
#                            stops with status 26, so an unattended series over an
#                            image that never decodes the class ends (default:
#                            unset, no limit)
#
# Exit status:
#   0  capture written    2  usage       3  environment or build
#   1  a self-test assertion did not hold
#   4  empty stream
#   5  no fingerprint line, or on the NUCLEO-H753ZI leg no boot handshake line
#      or one whose wait expired: either way the report's head may be lost
#   6  fingerprint mismatch (the resident image is not the built one)
#   7  no verdict line    8  verdict FAIL
#   9  no allocation-sensor canary line, or on the ESP32 leg no heap-view
#      canary line or one without its self-count
#  10  canary FAIL (the sensor did not observe its deliberate allocations), or
#      on the ESP32 leg a heap-view canary that is not PASS
#  11  a family's verdict line is missing or repeated, or on the ESP32 leg the
#      predictive family is not declared ABSENT on exactly one line
#  12  a family's verdict is not PASS
#  13  the image's kind line is missing or names the other image
#  14  no posture line, a posture field missing, or on the ESP32 leg no
#      schedule line
#  15  timing UNAVAILABLE: the cycle counter did not prove itself live, or the
#      timing image printed no timing line at all
#  16  a timing line without its sample count, or a per-iteration figure
#      without the block size it was divided by
#  17  a timed region lacks exactly one cold and one warm line
#  18  a stack watermark line reports a nonzero harness floor: the gap is
#      smaller than the harness's own frames, so its figure is not a
#      measurement of the chain
#  19  a stack watermark line reports a saturated window: the chain reached the
#      window's bottom word, so its depth is unknown rather than small (on the
#      ESP32 leg, a control-task high-water with zero bytes free)
#  20  a stack watermark line lacks its floor, its saturation verdict, its
#      figure or its reserve, or a measured region's line is missing or repeated
#      (on the ESP32 leg, the line's used and free do not sum to its stack)
#  21  a stack line's high-water exceeds the image's stack reserve: the
#      measured stack needs more than the linker reserved for it
#  22  the reserve a stack line prints disagrees with the _Min_Stack_Size
#      symbol in the ELF, so the line and the image disagree on the reserve (on
#      the ESP32 leg, with the control task's stack as the source names it)
#  23  one of the five stream tags is missing, so a block of the report is
#      absent rather than empty
#  24  a field name appears in both the exact-class and the spread-class
#      header, or twice in one, so it would be claimed bitwise identical and
#      bounded by a spread at once
#  25  a class block is malformed: its header does not key its rows (the loop
#      index for the exact class, the region for the spread class), a row's
#      width differs from its header's, an exact row's index breaks the
#      sequence 0, 1, 2, ... (a reordering is a difference, not a pass), or an
#      exact value is not a bit pattern of the leg's scalar: 16 hexadecimal
#      digits for binary64, 8 for binary32
#  26  a run's decoded reset class is not the class the series requires, or on
#      the ESP32 leg its header records a restart other than the one the
#      series' class names, so an EN-pin run is not pooled with power-on runs,
#      or a series kept COMPARE_MAX_REFUSED runs of another class
#  27  the exact-class blocks of two runs differ
#  28  a cycle field's spread across the runs exceeds the stated bound
#  29  the runs do not share what a comparison needs them to share: the build
#      id, the posture outside the reset fields, the spread rows and their
#      shape fields, or a spread column this script can classify
#  30  on the ESP32 leg, an armed window's heap view is not clean: a verdict
#      line carrying heap figures carries no heap_view verdict, or one other
#      than clean, so the window kept or dipped by more than the view's own
#      snapshot
#
# --- What this capture does NOT establish -------------------------------------
#
# Stated here rather than left implicit, because a gate that hides its blind
# spots is the same defect in a different place.
#
#   1. Copying the binary onto the mass-storage mount and calling sync confirms
#      a WRITE, not a PROGRAM. The mount's own bootloader does the programming
#      and reports nothing back over this path. The fingerprint check below is
#      what turns "a file was written" into "the target is running this image";
#      without it a successful copy proves nothing about the target.
#
#   2. The fingerprint is the linker's SHA-1 build id. It is an identity check
#      against ACCIDENTAL staleness -- a capture taken without re-flashing, a
#      board holding last week's image. It is NOT a signature and NOT
#      tamper-evidence: SHA-1 is collision-attackable, the note travels in the
#      clear inside the image, and nothing here authenticates who produced it.
#
#   3. A passing capture says the report arrived and matched. This script
#      cannot tell a power cycle from a warm reset; the image decodes its own
#      reset status register and prints the class on the posture line, and this
#      script requires and transports that label without judging it. The flags
#      accumulate until cleared and a power-on sets every flag a pin reset sets,
#      so the label describes this run only if the image that ran before it
#      cleared them: the first run after an image that never cleared the
#      register reports every reset since the last power-on, and reads as a
#      power-on.
#
#   4. The read is bounded, so a report that is slower than the timeout is
#      truncated. A truncated report fails on its missing verdict line rather
#      than passing on what did arrive, but the timeout is a limit chosen by this
#      script and not a property of the target.
#
#   5. A passing canary proves the allocation sensor observed one deliberate
#      allocation through each path it watches: a direct C allocation, a C++
#      new, and an Eigen heap vector. It does not extend that sight to the C
#      library's reentrant allocator entry points, which stdio uses and which
#      never pass through the wrapped symbols.
#
#   6. A spread row's figures are raw: each carries the instrument's own
#      overhead, which the report publishes on its own lines rather than
#      subtracting. The cold figure is one sample per boot; its spread is a
#      property of repeated captures, not of one.
#
#   7. A stack line's figure is the deepest word the chain disturbed in a
#      painted window, at word resolution and no finer than the gap on the same
#      line. A chain that writes a value equal to a word's pattern would hide
#      that word; the pattern carries each word's address, so no single value
#      can hide more than one. The compiler's per-function frames quoted in the
#      capture header attribute nothing to callees and are a lower bound. The
#      reserve refusal compares the reserve with the high-water of this run
#      only: a deeper path the run did not take is not covered, and nothing on
#      the target stops the stack at the reserve.
#
#   8. A family's PASS says the value its last step returned lies within its
#      derived bound of the host reference; it does not compare the steps before
#      it. The predictive family's value is the cost accumulated over the whole
#      run, so it does weigh every step, but its bound holds only under the
#      premise its own [verdict] predictive line reports (every solve stopped on the
#      stationarity test). The allocation count on the same line is
#      transported, not judged.
#
#   9. A series that compares identical shows the runs it holds were identical,
#      no more: zero differences in N runs leaves a per-run difference rate up
#      to about 3/N plausible at 95 percent confidence, so the run count is
#      part of the result. Every run restarts one resident image; a reflash
#      between runs is not part of the comparison, and only the [exact] block
#      is compared byte for byte, not the verdict lines.
#
#  10. The power-on prompt cannot verify the cycle it asks for. The node
#      vanishing proves the ST-LINK lost its USB connection, not that the
#      target's supply decayed; the decoded class is the only evidence, and it
#      is only as good as the priming run that cleared the flags before it.
#      The reader attaches after the ST-LINK enumerates; the image's boot
#      handshake holds the report until this script's byte arrives, and only
#      up to its timeout, so a replug slower than that still loses the head
#      and is refused with status 5. A series counts only runs that decoded
#      the requested class, so its exit status alone does not say how many
#      replugs it took: the refused runs are listed in its summary.
#      The handshake byte is a control character a modem manager's AT probe
#      never sends, but any process that writes that byte to the port during
#      the wait releases it as this script would.
#
#  11. The spread bound is one absolute cycle count applied to every timed
#      region's raw figures, from the 114-cycle overhead to the predictive
#      solve; it is the invocation's stated claim, not a figure derived here.
#
#  12. The self-test proves the parser, the fingerprint comparison, the class
#      assignment and every refusal against recordings and corruptions of them.
#      It proves nothing about the board: it opens no port, flashes nothing, and
#      its references are the values recorded beside the fixtures, not read from
#      a built image. The recordings are a fixed sample of the format: if the
#      firmware's report changes, the self-test goes on proving the parser
#      against the old format until a new recording replaces the fixture.
#
#  13. "discover" reports absence; it does not prove it. A board on a port the
#      discovery does not search reads as absent, and a SKIP line from a machine
#      that does have a board means the discovery, not the board, is missing.
#
#  On the ESP32 leg, in addition:
#
#  14. The fingerprint is the sha256 digest of the whole ELF, debug sections
#      included, so two builds of identical code from different paths carry
#      different ids; it is stricter than the code it names, never looser. Like
#      the NUCLEO's build id it is an identity check against staleness and not
#      a signature: this image enables no secure boot, and the descriptor that
#      carries the digest travels in the clear. esptool's verification after
#      writing proves the flash holds the bytes sent, which the mass-storage
#      copy of limit 1 cannot.
#
#  15. The reset label cannot tell a power cycle from an EN-pin reset: the part
#      latches the same chip power-on code for both, so the framework reports
#      ESP_RST_POWERON for a power cycle, the EN button, and the RTS-driven EN
#      pulse this script and esptool issue. Every capture this script restarts
#      therefore reads reset=power_on, and that label says nothing about the
#      supply. An en_pin series is named for the restart this script made, not
#      for the label, and is not power-on evidence. In --power-cycle mode the
#      console node vanishing proves the USB bridge lost its connection, not
#      that the supply decayed: a telemetry adapter wired to the board can feed
#      it through a signal pin. The labels this leg can separate are panic,
#      software, the watchdogs and brownout.
#
#  16. The stack figure is the control task's high-water as FreeRTOS counts it,
#      the untouched bytes of the fill pattern at the bottom of the task's
#      stack, read once after the report's last step. A value that equals the
#      fill byte hides itself, so the figure is a lower bound on the depth by
#      at most that run of bytes. It is judged against the task's stack as the
#      source constant names it, not against a symbol in the image, so a
#      source edited after the build would not be caught by this check; the
#      fingerprint is what ties the report to the build. The end-of-stack
#      watchpoint the posture reports costs up to 60 bytes of every task's
#      usable stack, so the figure is only comparable with one taken under the
#      same setting. The per-dimension stack probe prints on the console after
#      the report; its lines are transported in the side file and judged by
#      nothing here.
#
#  17. The allocation canary proves the counter saw a plain and an aligned C++
#      new; this leg does not wrap the C allocation family, so a C allocation
#      is visible only to the free-heap view printed beside it. That view's
#      own start call allocates a snapshot inside every window, so it never
#      reads zero: the image measures that self-count with an empty window,
#      prints it on the heap-view canary line, and judges each window against
#      it, and this script refuses a heap-view canary that did not see the
#      requested bytes beyond the self-count and a window that kept or dipped
#      by more than it. The view counts the whole heap, so another task's
#      allocation inside a window fails it too; it cannot say whose it was.
#
#  18. The timed task runs with the scheduler and interrupts enabled, so a
#      timing figure includes whatever preempted it. The interruption lines
#      count the samples the scheduler tick landed in, a lower bound on how
#      often a sample was interrupted; other interrupts leave no mark. The
#      timing image stops the task watchdog before it reports, since a timed
#      block outlasts its timeout without yielding, and its schedule line
#      prints the watchdog's state; this script transports that field and does
#      not require it.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

leg_dir="${repo_root}/examples/embedded/nucleo_h753zi"
build_dir="${CAPTURE_BUILD_DIR:-${repo_root}/build-arm-capture}"
mount_dir="${CAPTURE_MOUNT:-/run/media/${USER}/NOD_H753ZI}"
read_timeout="${CAPTURE_READ_TIMEOUT:-45}"
byte_cap="${CAPTURE_BYTE_CAP:-262144}"
build_jobs="${CAPTURE_BUILD_JOBS:-6}"

toolchain_prefix="${CTRLPP_ARM_TOOLCHAIN_PREFIX:-arm-none-eabi-}"
evidence_target="ctrlpp_nucleo_h753zi"
timing_target="ctrlpp_nucleo_h753zi_timing"
evidence_elf="${build_dir}/${evidence_target}.elf"

# Spelled out because ninja 1.13.2 aborts on this project's dynamic dependency
# file. Parallelism belongs to the build only; the capture is serial.
generator="Unix Makefiles"

# The image's closing line. Seeing it means the report is complete, so the read
# can stop before the timeout rather than always paying it.
report_end_marker="golden diff "

# Every tag the report carries, each owed at least one line.
expected_tags="meta canary exact spread verdict"

# Every family the image runs. Each must print exactly one verdict line, so a
# family that ran and reported nothing is a refusal rather than a shorter report.
expected_families="control estimation dsp trajectory predictive"

# Every region the evidence image measures the stack of, each owed one line.
expected_stack_regions="startup predictive.construct predictive.first_solve predictive.steady program"

# Every region the timing image times, each owed one cold and one warm line.
expected_regions="overhead control.design control.step estimation.step dsp.step trajectory.step predictive.step"

# Families the image declares absent, each owed exactly one ABSENT line.
absent_families=""

# The hexadecimal digits of one exact value: the scalar's bit width over four.
exact_digits=16

# Every posture field a report must carry, beside the reset class.
posture_fields="reset reset_flags reset_rsr core_hz icache dcache flash_latency_ws cycle_wrap_ms compiler"

# The reset fields a posture line carries, quoted per run and set aside when two
# runs' postures are compared.
reset_quote='reset=[^ ]* reset_flags=[^ ]* reset_rsr=[^ ]*'

# The reset classes "compare" can produce on this station.
collectable_classes="pin power_on"

# The NUCLEO-H753ZI image waits at boot for the host's byte; the ESP32 leg's
# telemetry adapter stays attached through a restart, so it needs none.
handshake_required="yes"

# The ESP32 evidence image prints a free-heap view beside its allocation count.
heap_view_judged="no"

board="nucleo"
board_label="NUCLEO-H753ZI"
console_tail="${CAPTURE_CONSOLE_TAIL:-15}"

# --- Arguments ----------------------------------------------------------------

subcommand=""
capture_mode="flash"
image_kind="evidence"
destination=""
replay="no"
self_test="no"
spread_bound=""
operands=()

usage()
{
    echo "  Usage: scripts/board_capture.sh [--board=esp32] prepare" >&2
    echo "         scripts/board_capture.sh [--board=esp32] capture [--timing] [--reset-only | --no-reset | --power-cycle] <destination-file>" >&2
    echo "         scripts/board_capture.sh [--board=esp32] check [--timing] <transcript-or-capture-file>" >&2
    echo "         scripts/board_capture.sh [--board=esp32] compare [--timing --bound=<cycles>] <count> <reset-class> <destination-directory>" >&2
    echo "         scripts/board_capture.sh [--board=esp32] compare --replay [--timing --bound=<cycles>] <reset-class> <series-directory>" >&2
    echo "         scripts/board_capture.sh [--board=esp32] discover" >&2
    echo "         scripts/board_capture.sh --self-test" >&2
}

# The ESP32 leg emits the same stream from a binary32 kernel without the
# predictive backend, on a part whose reset decode and stack figure differ.
select_board()
{
    board="$1"
    case "${board}" in
        nucleo)
            return 0
            ;;
        esp32)
            ;;
        *)
            echo "ERROR: unknown board '${board}'; expected 'esp32', or no --board for the NUCLEO-H753ZI." >&2
            exit 2
            ;;
    esac
    board_label="ESP32-POE-EA"
    leg_dir="${repo_root}/examples/embedded/esp32"
    build_dir="${CAPTURE_BUILD_DIR:-${leg_dir}/build}"
    evidence_target="ctrlpp_esp32"
    timing_target="ctrlpp_esp32_timing"
    evidence_elf="${build_dir}/${evidence_target}.elf"
    expected_families="control estimation dsp trajectory"
    absent_families="predictive"
    expected_stack_regions="program"
    expected_regions="overhead control.design control.step estimation.step dsp.step trajectory.step"
    exact_digits=8
    posture_fields="reset reset_reason reset_rom core_mhz cycle_bits cycle_wrap_ms stack_watchpoint assertions optimize compiler idf"
    reset_quote='reset=[^ ]* reset_reason=[^ ]* reset_rom=[^ ]*'
    collectable_classes="power_on en_pin"
    handshake_required="no"
    heap_view_judged="yes"
}

for arg in "$@"; do
    case "${arg}" in
        -h|--help)
            sed -n '2,/^$/p' "${BASH_SOURCE[0]}"
            exit 0
            ;;
        --timing)
            image_kind="timing"
            ;;
        --reset-only)
            capture_mode="reset-only"
            ;;
        --no-reset)
            capture_mode="no-reset"
            ;;
        --power-cycle)
            capture_mode="power-cycle"
            ;;
        --board=*)
            board="${arg#--board=}"
            ;;
        --replay)
            replay="yes"
            ;;
        --self-test)
            self_test="yes"
            ;;
        --bound=*)
            spread_bound="${arg#--bound=}"
            ;;
        -*)
            echo "ERROR: unknown option '${arg}'." >&2
            usage
            exit 2
            ;;
        *)
            if [ -z "${subcommand}" ]; then
                subcommand="${arg}"
            else
                operands+=("${arg}")
            fi
            ;;
    esac
done

if [ "${subcommand}" != "compare" ]; then
    if [ "${#operands[@]}" -gt 1 ]; then
        echo "ERROR: unexpected extra argument '${operands[1]}'." >&2
        exit 2
    fi
    destination="${operands[0]:-}"
fi

select_board "${board}"
case "${capture_mode}" in
    power-cycle)
        [ "${board}" = "esp32" ] \
            || { echo "ERROR: --power-cycle is an ESP32 capture mode; the NUCLEO-H753ZI's power-on runs are collected by 'compare'." >&2; exit 2; }
        ;;
esac

target_name="${evidence_target}"
if [ "${image_kind}" = "timing" ]; then
    target_name="${timing_target}"
fi
elf_path="${build_dir}/${target_name}.elf"
bin_path="${build_dir}/${target_name}.bin"

fail()
{
    local status="$1"
    shift
    echo "ERROR: $*" >&2
    exit "${status}"
}

# --- prepare ------------------------------------------------------------------

do_prepare()
{
    if [ "${board}" = "esp32" ]; then
        prepare_esp32
        return 0
    fi
    local configure_args=()
    read -r -a configure_args <<< "${CAPTURE_CONFIGURE_ARGS:-}"
    cmake -S "${leg_dir}" -B "${build_dir}" \
        -G "${generator}" \
        -DCMAKE_TOOLCHAIN_FILE="${leg_dir}/cmake/toolchain-arm-none-eabi.cmake" \
        -DCMAKE_BUILD_TYPE=Release \
        "${configure_args[@]}" \
        || fail 3 "configure of ${leg_dir} failed."

    cmake --build "${build_dir}" --parallel "${build_jobs}" \
        || fail 3 "build of ${leg_dir} failed."

    local target elf bin
    for target in "${evidence_target}" "${timing_target}"; do
        elf="${build_dir}/${target}.elf"
        bin="${build_dir}/${target}.bin"
        [ -f "${elf}" ] || fail 3 "the build produced no ${elf}."
        "${toolchain_prefix}objcopy" -O binary "${elf}" "${bin}" \
            || fail 3 "objcopy could not produce ${bin}."
        echo "prepared: ${elf}"
        echo "prepared: ${bin}"
        echo "build id: $(elf_build_id "${elf}")"
    done
}

# The framework owns this leg's build tree, its compiler and its image
# conversion, so the build is the framework's and only its outputs are checked.
prepare_esp32()
{
    [ -n "${IDF_PATH:-}" ] || fail 3 "the framework is not exported; run '. \$IDF_PATH/export.sh' (for example '. /opt/esp-idf/export.sh') first."
    [ -f "${build_dir}/CMakeCache.txt" ] \
        || fail 3 "${build_dir} is not configured; run 'idf.py -C ${leg_dir} -B ${build_dir} -G \"${generator}\" set-target esp32' once first."
    cmake --build "${build_dir}" --parallel "${build_jobs}" \
        || fail 3 "build of ${leg_dir} failed."

    local target elf bin
    for target in "${evidence_target}" "${timing_target}"; do
        elf="${build_dir}/${target}.elf"
        bin="${build_dir}/${target}.bin"
        [ -f "${elf}" ] || fail 3 "the build produced no ${elf}."
        [ -f "${bin}" ] || fail 3 "the build produced no ${bin}."
        echo "prepared: ${elf}"
        echo "prepared: ${bin}"
        echo "build id: $(elf_build_id "${elf}")"
        check_descriptor "${bin}" "$(elf_build_id "${elf}")"
    done
}

# The binary carries the digest of the ELF it was converted from; reading it
# back proves the binary about to be written is the ELF the report is checked
# against.
check_descriptor()
{
    local described=""
    described="$(python -m esptool --chip esp32 image-info "$1" 2>/dev/null \
        | sed -n 's/^ELF file SHA256:[[:space:]]*\([0-9a-f]\{64\}\)[[:space:]]*$/\1/p' | head -1)" \
        || fail 3 "esptool could not read the descriptor of $1; export the framework first."
    [ -n "${described}" ] || fail 3 "the descriptor of $1 names no ELF digest, or esptool is not available."
    [ "${described}" = "$2" ] \
        || fail 3 "$1 describes the ELF ${described}, not the built ${2}; the binary is stale against its ELF. Rebuild."
}

# --- discovery ----------------------------------------------------------------

# The first by-id link matching a pattern, resolved to its device node.
by_id_device()
{
    local pattern="$1" link
    for link in /dev/serial/by-id/${pattern}; do
        if [ -e "${link}" ]; then
            readlink -f "${link}"
            return 0
        fi
    done
    return 1
}

# The ESP32-POE-EA's own bridge is a CH340; CP210x bridges are the common
# alternative on other ESP32 boards. The telemetry adapter is an FT232R.
resolve_console_port()
{
    if [ -n "${CAPTURE_CONSOLE_PORT:-}" ]; then
        printf '%s\n' "${CAPTURE_CONSOLE_PORT}"
        return 0
    fi
    by_id_device '*1a86*' || by_id_device '*CP210*'
}

resolve_port()
{
    if [ -n "${CAPTURE_PORT:-}" ]; then
        printf '%s\n' "${CAPTURE_PORT}"
        return 0
    fi
    if [ "${board}" = "esp32" ]; then
        by_id_device '*FTDI*'
        return
    fi

    local link
    for link in /dev/serial/by-id/*STLINK*; do
        if [ -e "${link}" ]; then
            readlink -f "${link}"
            return 0
        fi
    done

    local device tty
    for device in /sys/bus/usb/devices/*/; do
        [ -r "${device}idVendor" ] || continue
        [ -r "${device}idProduct" ] || continue
        [ "$(cat "${device}idVendor")" = "0483" ] || continue
        [ "$(cat "${device}idProduct")" = "374e" ] || continue
        for tty in "${device}"*/tty/*; do
            if [ -e "${tty}" ]; then
                printf '/dev/%s\n' "$(basename "${tty}")"
                return 0
            fi
        done
    done

    return 1
}

# The linker script defines it as an absolute symbol, so nm prints its value.
# awk reads nm to the end: stopping early would kill nm with SIGPIPE, which
# pipefail turns into a failed assignment.
elf_stack_reserve()
{
    local value=""
    if [ "${board}" = "esp32" ]; then
        source_stack_reserve
        return 0
    fi
    value="$("${toolchain_prefix}nm" "${elf_path}" | awk '$3 == "_Min_Stack_Size" { print $1 }')"
    [ -n "${value}" ] || return 0
    printf '%d\n' "0x${value}"
}

# The control task's stack as the source names it. The expression is checked to
# be digits and multiplications before the shell evaluates it.
source_stack_reserve()
{
    local expression=""
    expression="$(sed -n 's/^inline constexpr std::uint32_t kControlTaskStackBytes = \([0-9][0-9 *]*\);$/\1/p' "${leg_dir}/main/stack_probe.h")"
    printf '%s\n' "${expression}" | grep -Eqx '[0-9]+( \* [0-9]+)*' || return 0
    echo "$((expression))"
}

# The compiler the build used, as the framework recorded it.
toolchain_line()
{
    local compiler=""
    if [ "${board}" != "esp32" ]; then
        "${toolchain_prefix}gcc" --version | head -1
        return 0
    fi
    compiler="$(sed -n 's/^ *"c_compiler": *"\([^"]*\)".*/\1/p' "${build_dir}/project_description.json" 2>/dev/null | head -1)"
    [ -x "${compiler}" ] || { echo "unknown (no c_compiler in ${build_dir}/project_description.json)"; return 0; }
    "${compiler}" --version | head -1
}

# The framework revision the build was configured against.
framework_line()
{
    sed -n 's/^ *"git_revision": *"\([^"]*\)".*/ESP-IDF \1/p' "${build_dir}/project_description.json" 2>/dev/null | head -1
}

# What a transcript is checked against: the fingerprint, the stack reserve and
# the toolchain of the image it claims to come from. Read from the built ELF
# here; the self-test reads the values recorded beside its fixtures instead.
reference_id=""
reference_reserve=""
reference_origin=""
reference_evidence_id=""
reference_toolchain=""
reference_predates=""

# A recording made before the image printed a line names that line in its
# reference, so the self-test replays it under the checks it was made for. A
# built image predates nothing.
predates()
{
    case " ${reference_predates} " in
        *" $1 "*) return 0 ;;
    esac
    return 1
}

load_elf_reference()
{
    [ -f "${elf_path}" ] || fail 3 "no ${elf_path}; run 'board_capture.sh prepare' first."
    reference_id="$(elf_build_id)"
    [ -n "${reference_id}" ] || fail 3 "${elf_path} carries no NT_GNU_BUILD_ID note; the image cannot be fingerprinted."
    reference_predates="$(elf_predates)"
    reference_reserve="$(elf_stack_reserve)"
    reference_origin="${elf_path}"
    reference_toolchain="$(toolchain_line)"
    if [ "${image_kind}" = "timing" ]; then
        reference_evidence_id="$(elf_build_id "${evidence_elf}")"
    fi
}

# The report lines an image predates, read from its symbols: an image built
# before a line existed carries no code to print it, so its transcripts are
# checked as they were made. Here-strings rather than pipes, so grep stopping at
# its first match cannot fail the pipeline under pipefail.
elf_predates()
{
    local nm="${toolchain_prefix}nm" compiler="" symbols=""
    if [ "${board}" = "esp32" ]; then
        compiler="$(sed -n 's/^ *"c_compiler": *"\([^"]*\)".*/\1/p' "${build_dir}/project_description.json" 2>/dev/null | head -1)"
        nm="${compiler%gcc}nm"
    fi
    symbols="$("${nm}" "${elf_path}" 2>/dev/null)" || fail 3 "'${nm}' could not read the symbols of ${elf_path}, so the report lines the image prints are unknown."
    if [ "${handshake_required}" = "yes" ] && ! grep -q 'await_boot_handshake' <<< "${symbols}"; then
        printf 'boot_handshake '
    fi
    if [ "${heap_view_judged}" = "yes" ] && ! grep -q 'alloc_sensor_heap_self' <<< "${symbols}"; then
        printf 'heap_self '
    fi
}

elf_build_id()
{
    if [ "${board}" = "esp32" ]; then
        sha256sum "${1:-${elf_path}}" | cut -d' ' -f1
        return 0
    fi
    "${toolchain_prefix}readelf" -n "${1:-${elf_path}}" \
        | sed -n 's/^[[:space:]]*Build ID:[[:space:]]*\([0-9a-f]\{40\}\)[[:space:]]*$/\1/p' \
        | head -1
}

# --- the reader ---------------------------------------------------------------
#
# cat rather than head -c: cat writes through without stdio buffering, so the
# transcript is greppable while the report is still arriving and the read can
# stop at the verdict line. The byte cap is enforced by the polling loop below.

reader_pid=""
raw_stream=""
transcript=""

cleanup()
{
    local scratch=""
    if [ -n "${reader_pid}" ]; then
        kill "${reader_pid}" 2>/dev/null || true
        wait "${reader_pid}" 2>/dev/null || true
        reader_pid=""
    fi
    console_stop
    keep_raw
    for scratch in "${console_raw}" "${flash_log}"; do
        [ -z "${scratch}" ] || rm -f -- "${scratch}"
    done
    if [ -n "${raw_stream}" ]; then
        rm -f -- "${raw_stream}"
    fi
    if [ -n "${transcript}" ]; then
        rm -f -- "${transcript}"
    fi
    return 0
}

keep_raw()
{
    local stamp="" scratch=""
    [ -n "${CAPTURE_KEEP_RAW:-}" ] || return 0
    mkdir -p -- "${CAPTURE_KEEP_RAW}" || return 0
    stamp="$(date -u '+%Y%m%dT%H%M%S')-$$"
    for scratch in "${raw_stream}" "${console_raw}"; do
        [ -z "${scratch}" ] || [ ! -f "${scratch}" ] \
            || cp -- "${scratch}" "${CAPTURE_KEEP_RAW}/${stamp}-$(basename "${scratch}")"
    done
    return 0
}

reader_start()
{
    local port="$1"
    stty -F "${port}" 115200 raw -echo -echoe -echok
    timeout "${read_timeout}" cat < "${port}" >> "${raw_stream}" &
    reader_pid=$!
}

reader_stop()
{
    if [ -n "${reader_pid}" ]; then
        kill "${reader_pid}" 2>/dev/null || true
        wait "${reader_pid}" 2>/dev/null || true
        reader_pid=""
    fi
}

# The mass-storage copy momentarily re-enumerates the ST-LINK CDC, which tears
# down a live reader's descriptor. Wait for the node to go and come back, then
# re-attach; the target boots and streams a few seconds after the copy, so the
# re-attached reader still catches the whole boot-once report.
wait_for_port_cycle()
{
    local port="$1" waited=0
    while [ "${waited}" -lt 10 ] && [ -e "${port}" ]; do
        sleep 1
        waited=$((waited + 1))
    done
    waited=0
    while [ "${waited}" -lt "${read_timeout}" ] && [ ! -e "${port}" ]; do
        sleep 1
        waited=$((waited + 1))
    done
    [ -e "${port}" ] || fail 3 "the serial port ${port} did not come back after programming."
    sleep 1
}

# The debug interface intermittently refuses the first connect while the console
# is still settling after a re-enumeration, reporting a core_id read failure. A
# single attempt therefore fails a capture for a reason that has nothing to do
# with the board's report, so the connect is retried a bounded number of times.
board_reset()
{
    local attempt=1
    while [ "${attempt}" -le 3 ]; do
        if st-flash reset; then
            return 0
        fi
        sleep 2
        attempt=$((attempt + 1))
    done
    fail 3 "st-flash reset did not connect to the target in 3 attempts; it was not restarted."
}

wait_for_report()
{
    local port="$1" waited=0 size=0
    while [ "${waited}" -lt "${read_timeout}" ]; do
        size="$(wc -c < "${raw_stream}" 2>/dev/null || echo 0)"
        if [ "${size}" -ge "${byte_cap}" ]; then
            return 0
        fi
        if grep -q "${report_end_marker}" "${raw_stream}" 2>/dev/null; then
            return 0
        fi
        release_handshake "${port}"
        sleep 1
        waited=$((waited + 1))
    done
    return 0
}

# The byte the image waits for, as boot_handshake.h names it, checked to be one
# hexadecimal byte and turned into the octal escape POSIX printf takes.
handshake_escape=""

load_handshake_byte()
{
    local value=""
    value="$(sed -n 's/^inline constexpr std::uint8_t kBootHandshakeByte = 0x\([0-9A-Fa-f][0-9A-Fa-f]\);$/\1/p' "${leg_dir}/boot_handshake.h")"
    [ -n "${value}" ] || fail 3 "${leg_dir}/boot_handshake.h names no kBootHandshakeByte, so the byte the image waits for at boot is unknown."
    printf -v handshake_escape '\\%03o' "0x${value}"
}

# Sent until the report opens: a byte sent while the target is still held in
# reset is lost, and once the image has taken one it no longer reads the line.
# The open is bounded because a serial node can block an open for its carrier.
release_handshake()
{
    [ "${handshake_required}" = "yes" ] || return 0
    if grep -q '^\[meta\] build_id=' "${raw_stream}" 2>/dev/null; then
        return 0
    fi
    timeout 2 sh -c 'printf "$1" > "$2"' sh "${handshake_escape}" "$1" 2>/dev/null || true
}

# --- refusals -----------------------------------------------------------------

verify_transcript()
{
    local printed_id=""

    grep -q '[^[:space:]]' "${transcript}" \
        || fail 4 "the capture read an EMPTY stream from the board -- no report was observed, and an absent report is not a pass. The image reports once at boot: reset the target while this script holds the port (drop the --no-reset flag), or press B2 and re-run."

    printed_id="$(sed -n 's/^\[meta\] build_id=\([0-9a-f]*\).*$/\1/p' "${transcript}" | head -1)"
    [ -n "${printed_id}" ] \
        || fail 5 "the captured report carries no '[meta] build_id=' line, so the image that produced it cannot be identified."

    [ "${printed_id}" = "${reference_id}" ] \
        || fail 6 "the running image is NOT the built one: the board reported ${printed_id} but ${reference_origin} is ${reference_id}. Flash the target before capturing."

    verify_handshake
    verify_kind
    verify_posture

    if [ "${image_kind}" = "evidence" ]; then
        verify_canary
    fi

    grep -q "${report_end_marker}" "${transcript}" \
        || fail 7 "the captured report carries no '${report_end_marker}' verdict line -- it is truncated or the run did not complete."

    if grep -q "${report_end_marker}FAIL" "${transcript}"; then
        fail 8 "the board's golden verdict is FAIL."
    fi

    verify_families
    if [ "${image_kind}" = "evidence" ] && [ "${heap_view_judged}" = "yes" ]; then
        verify_heap_windows
    fi
    verify_tags
    verify_spread_header

    if [ "${image_kind}" = "timing" ]; then
        verify_timing
    fi
    if [ "${image_kind}" = "evidence" ] || [ "${board}" = "esp32" ]; then
        verify_stack
    fi

    verify_classes
    verify_exact
    return 0
}

# The lines carrying one tag, with the tag removed.
block_lines()
{
    sed -n "s/^\[$1\] //p" "${transcript}"
}

# A report whose wait expired began with no reader attached, so its head may be
# missing even where the fingerprint line survived.
verify_handshake()
{
    local line=""
    [ "${handshake_required}" = "yes" ] || return 0
    if predates boot_handshake; then
        return 0
    fi
    line="$(grep -m 1 '^\[meta\] boot_handshake ' "${transcript}" || true)"
    [ -n "${line}" ] \
        || fail 5 "the captured report carries no '[meta] boot_handshake' line, so whether its head reached the reader is unstated."
    printf '%s\n' "${line}" | grep -q '^\[meta\] boot_handshake state=received ' \
        || fail 5 "the image's boot handshake was not received: '${line}'. It began reporting with no reader attached, so the head of its report may be lost."
    return 0
}

# The heap view never reads zero, since its own snapshot lands in every window,
# so its canary is live only if the window kept the requested bytes beyond that
# self-count, and the self-count must be printed for the windows to be judged.
verify_heap_canary()
{
    local line=""
    line="$(grep -m 1 '^\[canary\] heap_view ' "${transcript}" || true)"
    [ -n "${line}" ] \
        || fail 9 "the captured report carries no '[canary] heap_view' line, so the heap view printed beside the allocation figures never proved it can see an allocation."
    printf '%s\n' "${line}" | grep -q '^\[canary\] heap_view PASS ' \
        || fail 10 "the heap-view canary did not pass: '${line}'. A view that did not see the canary's bytes beyond its own snapshot is blind to the allocations it is printed beside."
    if predates heap_self; then
        return 0
    fi
    printf '%s\n' "${line}" | grep -Eq ' self=[0-9]+( |$)' \
        || fail 9 "the heap-view canary line carries no self-count, so the windows' heap figures have nothing to be judged against: '${line}'."
    return 0
}

# Each window's view is judged on the board against the self-count; one that
# kept or dipped by more saw an allocation, the library's or another task's.
verify_heap_windows()
{
    local line=""
    while IFS= read -r line; do
        if printf '%s\n' "${line}" | grep -q ' heap_view='; then
            printf '%s\n' "${line}" | grep -Eq ' heap_view=clean( |$)' \
                || fail 30 "an armed window's heap view is not clean: it kept or dipped by more than the view's own snapshot: '${line}'."
        elif ! predates heap_self; then
            fail 30 "a line carrying heap figures carries no heap_view verdict: '${line}'."
        fi
    done < <(grep '^\[verdict\] .* heap_retained=' "${transcript}" || true)
    return 0
}

verify_tags()
{
    local tag=""
    for tag in ${expected_tags}; do
        grep -q "^\[${tag}\] " "${transcript}" \
            || fail 23 "the captured report carries no '[${tag}]' line, so that block of the report is absent."
    done
    return 0
}

# A field in both headers would be asserted bitwise identical and bounded by a
# spread at once, which is the mixing the two classes exist to prevent.
verify_classes()
{
    local repeated=""
    repeated="$( { block_lines exact | head -1; block_lines spread | head -1; } | tr ',' '\n' | sort | uniq -d | head -1)"
    [ -z "${repeated}" ] \
        || fail 24 "the field '${repeated}' is named in both class headers, or twice in one: '$(block_lines exact | head -1)' and '$(block_lines spread | head -1)'."
    return 0
}

# Positional, never sorted: rows that agree as a set but not in order are a
# difference in the run, and a value that is not a bit pattern cannot carry a
# bitwise claim. The width is tested with length() rather than an interval
# expression, which older mawk builds do not parse.
verify_exact()
{
    local problem=""
    problem="$(block_lines exact | awk -F, -v digits="${exact_digits}" '
        NR == 1 { width = NF; if ($1 != "k") { print "the exact-class header does not key its rows by the loop index: \047" $0 "\047"; exit } next }
        NF != width { print "an exact-class row carries " NF " fields against a header of " width ": \047" $0 "\047"; exit }
        $1 != (NR - 2) "" { print "an exact-class row carries index " $1 " where " NR - 2 " is required: \047" $0 "\047"; exit }
        { for (i = 2; i <= NF; ++i) if (length($i) != digits || $i ~ /[^0-9A-F]/) { print "an exact-class value is not a bit pattern of " digits " hexadecimal digits: \047" $0 "\047"; exit } }
        END { if (NR == 1) print "the exact-class block carries a header and no rows." }')"
    [ -z "${problem}" ] || fail 25 "${problem}"
    return 0
}

verify_spread_header()
{
    [ "$(block_lines spread | head -1 | cut -d, -f1)" = "region" ] \
        || fail 25 "the spread-class header does not key its rows by region: '$(block_lines spread | head -1)'."
    return 0
}

# The header that keys the spread rows, which the timing refusals index by name.
spread_column()
{
    block_lines spread | head -1 | tr ',' '\n' | grep -n -x -m 1 "$1" | cut -d: -f1
}

# The two images run the same report, so a capture of the wrong one could
# otherwise pass every refusal the other image's mode applies.
verify_kind()
{
    local kind=""
    kind="$(sed -n 's/^\[meta\] kind=\([a-z]*\).*$/\1/p' "${transcript}" | head -1)"
    [ "${kind}" = "${image_kind}" ] \
        || fail 13 "the report names its image kind as '${kind:-absent}' where this capture requires '${image_kind}'."
    return 0
}

verify_posture()
{
    local posture="" field=""
    posture="$(grep -m 1 '^\[meta\] posture ' "${transcript}" || true)"
    [ -n "${posture}" ] \
        || fail 14 "the captured report carries no '[meta] posture' line, so the reset and the configuration its figures were taken under are unstated."
    for field in ${posture_fields}; do
        printf '%s\n' "${posture}" | grep -Eq " ${field}=[^ ]+" \
            || fail 14 "the posture line carries no '${field}=' field: '${posture}'."
    done
    if [ "${board}" = "esp32" ]; then
        grep -q '^\[meta\] schedule ' "${transcript}" \
            || fail 14 "the captured report carries no '[meta] schedule' line, so the scheduler, interrupt and code-placement posture its timing ran under is unstated."
    fi
    return 0
}

# A figure from a counter that never proved it runs is not a measurement, and a
# worst-of-N without its N, or a per-iteration figure without its divisor, is
# not a bound.
verify_timing()
{
    local line="" region="" regime="" count=0 rows="" width=0 n_at="" block_at=""
    line="$(grep -m 1 '^\[canary\] cycle_counter ' "${transcript}" || true)"
    printf '%s\n' "${line}" | grep -q '^\[canary\] cycle_counter PASS' \
        || fail 15 "the cycle counter did not prove itself live: '${line:-no [canary] cycle_counter line}'."
    rows="$(block_lines spread | tail -n +2)"
    [ -n "${rows}" ] || fail 15 "the timing image printed no spread-class row."

    width="$(block_lines spread | head -1 | awk -F, '{ print NF }')"
    n_at="$(spread_column n)"
    block_at="$(spread_column block)"
    while IFS= read -r line; do
        [ "$(printf '%s\n' "${line}" | awk -F, '{ print NF }')" -eq "${width}" ] \
            || fail 25 "a spread-class row's width differs from its header's: '${line}'."
        [ -n "${n_at}" ] && printf '%s\n' "${line}" | cut -d, -f"${n_at}" | grep -Eqx '[1-9][0-9]*' \
            || fail 16 "a timing row carries no sample count: '${line}'."
        [ -n "${block_at}" ] && printf '%s\n' "${line}" | cut -d, -f"${block_at}" | grep -Eqx '[1-9][0-9]*' \
            || fail 16 "a per-iteration figure carries no block size: '${line}'."
    done <<< "${rows}"

    for region in ${expected_regions}; do
        for regime in cold warm; do
            count="$(printf '%s\n' "${rows}" | grep -c "^${region},${regime}," || true)"
            [ "${count}" -eq 1 ] \
                || fail 17 "the report carries ${count} '${region}' ${regime} rows where exactly one is required."
        done
    done
    return 0
}

# A saturated window is not a smaller measurement but no measurement, a figure
# beside a nonzero floor measured the harness rather than the chain, and a
# high-water past the reserve ran on memory the link never set aside for it.
verify_stack()
{
    local line="" region="" count=0 reserve="" printed="" from_top=""
    if [ "${board}" = "esp32" ]; then
        verify_task_stack
        return 0
    fi
    reserve="${reference_reserve}"
    [ -n "${reserve}" ] || fail 3 "${reference_origin} names no stack reserve (_Min_Stack_Size), so the reserve is unknown."
    while IFS= read -r line; do
        printf '%s\n' "${line}" | grep -Eq ' floor=[0-9]+( |$)' \
            || fail 20 "a stack line carries no floor: '${line}'."
        printf '%s\n' "${line}" | grep -Eq ' saturated=(yes|no)( |$)' \
            || fail 20 "a stack line carries no saturation verdict: '${line}'."
        printf '%s\n' "${line}" | grep -Eq ' floor=0( |$)' \
            || fail 18 "a stack line reports a nonzero harness floor, so its figure is withheld: '${line}'."
        printf '%s\n' "${line}" | grep -q ' saturated=no' \
            || fail 19 "a stack line reports a saturated window, which is no measurement: '${line}'."
        printf '%s\n' "${line}" | grep -Eq ' used=[0-9]+ from_top=[0-9]+ free=[0-9]+ of=[0-9]+ ' \
            || fail 20 "a stack line with a zero floor and an unsaturated window carries no figure: '${line}'."
        printed="$(printf '%s\n' "${line}" | sed -n 's/.* reserve=\([0-9][0-9]*\) .*/\1/p')"
        [ -n "${printed}" ] \
            || fail 20 "a stack line carries no reserve: '${line}'."
        [ "${printed}" -eq "${reserve}" ] \
            || fail 22 "a stack line prints a reserve of ${printed} B where ${reference_origin} has _Min_Stack_Size = ${reserve} B: '${line}'."
        from_top="$(printf '%s\n' "${line}" | sed -n 's/.* from_top=\([0-9][0-9]*\) .*/\1/p')"
        [ "${from_top}" -le "${reserve}" ] \
            || fail 21 "the measured stack high-water of ${from_top} B exceeds the image's ${reserve} B reserve: '${line}'."
    done < <(grep '^\[verdict\] stack ' "${transcript}")

    for region in ${expected_stack_regions}; do
        count="$(grep -c "^\[verdict\] stack region=${region} " "${transcript}" || true)"
        [ "${count}" -eq 1 ] \
            || fail 20 "the report carries ${count} '${region}' stack lines where exactly one is required."
    done
    return 0
}

# A FreeRTOS high-water is the fill pattern left untouched at the bottom of the
# task's stack. Zero bytes left means the task reached its last byte, which is a
# depth this figure cannot state, and a figure whose parts do not sum to the
# stack it names is not a figure of that stack.
verify_task_stack()
{
    local line="" used="" free="" of=""
    [ -n "${reference_reserve}" ] || fail 3 "${leg_dir}/main/stack_probe.h names no control-task stack (kControlTaskStackBytes), so the stack is unknown."
    [ "$(grep -c '^\[verdict\] stack region=program ' "${transcript}" || true)" -eq 1 ] \
        || fail 20 "the report carries $(grep -c '^\[verdict\] stack region=program ' "${transcript}" || true) 'program' stack lines where exactly one is required."
    line="$(grep -m 1 '^\[verdict\] stack region=program ' "${transcript}")"
    printf '%s\n' "${line}" | grep -Eq ' used=[0-9]+ free=[0-9]+ of=[0-9]+ bytes$' \
        || fail 20 "the control task's stack line carries no figure: '${line}'."
    used="$(printf '%s\n' "${line}" | sed -n 's/.* used=\([0-9]*\) .*/\1/p')"
    free="$(printf '%s\n' "${line}" | sed -n 's/.* free=\([0-9]*\) .*/\1/p')"
    of="$(printf '%s\n' "${line}" | sed -n 's/.* of=\([0-9]*\) .*/\1/p')"
    [ "$((used + free))" -eq "${of}" ] \
        || fail 20 "the control task's stack line does not add up: used ${used} and free ${free} against ${of} bytes: '${line}'."
    [ "${of}" -eq "${reference_reserve}" ] \
        || fail 22 "the control task's stack line names a ${of} B stack where the source names ${reference_reserve} B: '${line}'."
    [ "${free}" -gt 0 ] \
        || fail 19 "the control task's stack was used to its last byte, so its depth is unknown rather than ${used} B: '${line}'."
    return 0
}

verify_families()
{
    local family="" line="" count=0
    for family in ${expected_families}; do
        count="$(grep -c "^\[verdict\] family ${family} " "${transcript}" || true)"
        [ "${count}" -eq 1 ] \
            || fail 11 "the captured report carries ${count} '[verdict] family ${family}' lines where exactly one is required."
        line="$(grep -m 1 "^\[verdict\] family ${family} " "${transcript}")"
        printf '%s\n' "${line}" | grep -q ' verdict=PASS ' \
            || fail 12 "the ${family} family did not pass: '${line}'."
    done
    for family in ${absent_families}; do
        count="$(grep -c "^\[verdict\] family ${family} " "${transcript}" || true)"
        line="$(grep -m 1 "^\[verdict\] family ${family} " "${transcript}" || true)"
        [ "${count}" -eq 1 ] && printf '%s\n' "${line}" | grep -Eq ' verdict=ABSENT( |$)' \
            || fail 11 "the ${family} family, which this image is built without, must be declared ABSENT on exactly one line; the report carries ${count}: '${line}'."
    done
    return 0
}

# The allocation figure is only evidence behind a canary that proved the sensor
# can see an allocation, so an absent or failing canary fails the whole capture.
verify_canary()
{
    local canary_line=""
    canary_line="$(grep -m 1 '^\[canary\] allocation ' "${transcript}" || true)"
    [ -n "${canary_line}" ] \
        || fail 9 "the captured report carries no '[canary] allocation' line, so its allocation figure comes from a sensor that never proved it can see an allocation."

    printf '%s\n' "${canary_line}" | grep -q '^\[canary\] allocation PASS observed=[1-9]' \
        || fail 10 "the allocation-sensor canary did not pass: '${canary_line}'. A zero reported by a sensor that missed its own deliberate allocations is blindness, not evidence."

    if [ "${heap_view_judged}" = "yes" ]; then
        verify_heap_canary
    fi
    return 0
}

write_artifact()
{
    local port="$1" built_id="$2" mount_note="not used in this mode"
    if [ "${capture_mode}" = "flash" ]; then
        mount_note="${mount_dir}"
    fi
    {
        echo "# ctrlpp NUCLEO-H753ZI on-target capture"
        echo "# captured:  $(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo "# mode:      ${capture_mode}"
        echo "# port:      ${port}"
        echo "# mount:     ${mount_note}"
        echo "# image:     ${elf_path}"
        echo "# kind:      ${image_kind}"
        echo "# build id:  ${built_id}"
        if [ "${image_kind}" = "timing" ]; then
            echo "# evidence:  $(elf_build_id "${evidence_elf}") (the allocation-evidence image built from the same sources)"
        fi
        echo "# toolchain: $("${toolchain_prefix}gcc" --version | head -1)"
        echo "#"
        echo "# The build id above is the linker's SHA-1 note, read out of the ELF with"
        echo "# ${toolchain_prefix}readelf -n and matched against the line the firmware printed."
        echo "# It is an identity check against accidental staleness, not a signature."
        echo "#"
        echo "# Handshake: $(sed -n 's/^\[meta\] boot_handshake //p' "${transcript}" | head -1), the image's"
        echo "# wait at boot for this script's byte before it printed anything."
        echo "#"
        echo "# Reset: $(sed -n 's/^\[meta\] posture \(reset=[^ ]* reset_flags=[^ ]* reset_rsr=[^ ]*\).*/\1/p' "${transcript}" | head -1),"
        echo "# as the image decoded its reset status register before clearing it. The label"
        echo "# describes this run only if the image that ran before it cleared the flags."
        echo "#"
        echo "# The [exact] rows print each double as its binary64 bit pattern in 16"
        echo "# uppercase hexadecimal digits, so two runs that differ in the last bit differ"
        echo "# in the row."
        if [ "${image_kind}" = "timing" ]; then
            band_cause
        else
            static_frames
            echo "#"
            echo "# Stack reserve: $(elf_stack_reserve) bytes, the _Min_Stack_Size symbol read"
            echo "# with ${toolchain_prefix}nm out of the image above. Every [verdict] stack line below"
            echo "# prints the same figure and none of their high-waters exceeds it."
        fi
        echo
        cat "${transcript}"
    } > "${destination}"
}

# The compiler's prediction, quoted beside the runtime figure it does not replace.
static_frames()
{
    local su_dir="${build_dir}/CMakeFiles/${target_name}.dir" worst=""
    worst="$(find "${su_dir}" -name '*.su' -exec cat {} + 2>/dev/null \
        | awk -F '\t' '$2 + 0 > worst { worst = $2 + 0; line = $0 } END { print line }')"
    echo "#"
    if [ -z "${worst}" ]; then
        echo "# Static frames: no -fstack-usage report was found under ${su_dir}."
        return 0
    fi
    echo "# Static frames: the worst named frame in this image's -fstack-usage report is"
    echo "# $(printf '%s\n' "${worst}" | cut -f2) bytes ($(printf '%s\n' "${worst}" | cut -f3)), in"
    echo "# $(printf '%s\n' "${worst}" | cut -f1 | sed -e "s|^${repo_root}/||" -e 's/ \[with .*//')."
    echo "# That report attributes nothing to callees, so a frame or a sum of frames is a"
    echo "# lower bound on what a stack must hold. The [verdict] stack lines below are the runtime"
    echo "# watermark over the whole chain, the figure a stack must actually cover. Both"
    echo "# are published; neither replaces the other."
}

# Assembled from the posture and work lines, so it can name a cache as a cause
# only when the image reports that cache enabled, and quotes the solver's own
# iteration counts rather than asserting what varied.
band_cause()
{
    local posture="" latency="" icache="" dcache="" caches="" cold="" warm=""
    posture="$(grep -m 1 '^\[meta\] posture ' "${transcript}")"
    latency="$(printf '%s\n' "${posture}" | sed -n 's/.* flash_latency_ws=\([^ ]*\).*/\1/p')"
    icache="$(printf '%s\n' "${posture}" | sed -n 's/.* icache=\([^ ]*\).*/\1/p')"
    dcache="$(printf '%s\n' "${posture}" | sed -n 's/.* dcache=\([^ ]*\).*/\1/p')"
    caches="both caches are off (icache=${icache} dcache=${dcache}), so no part of the cold-to-warm difference is a cache effect"
    if [ "${icache}" != "off" ] || [ "${dcache}" != "off" ]; then
        caches="the posture reports icache=${icache} dcache=${dcache}, so the cold-to-warm difference includes cache fill"
    fi
    cold="$(sed -n 's/^\[verdict\] work region=predictive.step regime=cold iterations=\([0-9]*\).*/\1/p' "${transcript}")"
    warm="$(sed -n 's/^\[verdict\] work region=predictive.step regime=warm .*iterations_min=\([0-9]*\) iterations_max=\([0-9]*\).*/\1 to \2/p' "${transcript}")"
    echo "#"
    echo "# Band cause: ${caches}. Beyond the caches it can come only from branch"
    echo "# prediction and instruction-fetch buffering in front of flash read at ${latency}"
    echo "# wait states, and this capture does not separate the two. The predictive"
    echo "# step's band is solver work: its cold solve ran ${cold:-?} iterations against"
    echo "# ${warm:-?} in the window."
}

# --- the ESP32 transport ------------------------------------------------------

console_raw=""
console_pid=""
console_file=""
flash_log=""
chip_line=""
stream_mark=0
console_mark=0
discarded=0

console_stop()
{
    if [ -n "${console_pid}" ]; then
        kill "${console_pid}" 2>/dev/null || true
        wait "${console_pid}" 2>/dev/null || true
        console_pid=""
    fi
}

# esptool leaves the port's termios with a zero minimum read, under which a
# blocking reader sees end of file at the first pause, so the reader restates it.
esp32_reader()
{
    stty -F "$1" 115200 raw -echo -echoe -echok min 1 time 0
    timeout "$3" cat < "$1" >> "$2" &
}

# The by-id link of the board's own USB bridge. The link survives a replug
# where the tty number may not, so a power cycle is watched through it.
console_link()
{
    local link
    if [ -n "${CAPTURE_CONSOLE_PORT:-}" ]; then
        printf '%s\n' "${CAPTURE_CONSOLE_PORT}"
        return 0
    fi
    for link in /dev/serial/by-id/*1a86* /dev/serial/by-id/*CP210*; do
        if [ -e "${link}" ]; then
            printf '%s\n' "${link}"
            return 0
        fi
    done
    return 1
}

# The lines are set through a descriptor of their own with ioctls, leaving the
# termios the readers rely on untouched. The kernel raises DTR and RTS on every
# open, which this board's reset circuit reads as run; DTR low with RTS high
# pulls EN low. The console reader holds the port open throughout, so closing
# this descriptor does not drop the lines.
drive_en()
{
    python3 - "$1" "$2" <<'PYTHON' || fail 3 "could not drive the reset lines of $1."
import fcntl, os, struct, sys, termios
fd = os.open(sys.argv[1], os.O_RDWR | os.O_NOCTTY | os.O_NONBLOCK)
if sys.argv[2] == "hold":
    fcntl.ioctl(fd, termios.TIOCMBIC, struct.pack("I", termios.TIOCM_DTR))
    fcntl.ioctl(fd, termios.TIOCMBIS, struct.pack("I", termios.TIOCM_RTS))
else:
    fcntl.ioctl(fd, termios.TIOCMBIS, struct.pack("I", termios.TIOCM_DTR | termios.TIOCM_RTS))
os.close(fd)
PYTHON
}

# The build's own flash_args, with the image of the requested kind in the
# application slot. The chip is left in the flasher, not restarted, so nothing
# runs until this script releases EN with its readers attached.
flash_esp32()
{
    local console="$1" settings="" offset="" file="" app="no" regions=0 args=()
    [ -f "${build_dir}/flash_args" ] || fail 3 "no ${build_dir}/flash_args; the framework's build writes it."
    settings="$(head -1 "${build_dir}/flash_args")"
    printf '%s\n' "${settings}" | grep -Eqx '(--flash-[a-z]+ [0-9A-Za-z]+ ?)+' \
        || fail 3 "the first line of ${build_dir}/flash_args is not a list of flash settings: '${settings}'."
    while read -r offset file; do
        [ -n "${offset}" ] || continue
        if [ "${file}" = "${evidence_target}.bin" ]; then
            file="$(basename "${bin_path}")"
            app="yes"
        fi
        args+=("${offset}" "${build_dir}/${file}")
        regions=$((regions + 1))
    done < <(tail -n +2 "${build_dir}/flash_args")
    [ "${app}" = "yes" ] || fail 3 "${build_dir}/flash_args names no ${evidence_target}.bin application slot."
    # shellcheck disable=SC2086
    python -m esptool --chip esp32 -p "${console}" -b 460800 --before default-reset --after no-reset \
        write-flash ${settings} "${args[@]}" > "${flash_log}" 2>&1 \
        || { cat "${flash_log}" >&2; fail 3 "esptool could not write ${bin_path} to the board on ${console}."; }
    [ "$(grep -c 'Hash of data verified' "${flash_log}" || true)" -eq "${regions}" ] \
        || { cat "${flash_log}" >&2; fail 3 "esptool did not report a verified hash for each of the ${regions} regions it wrote."; }
    chip_line="$(sed -n 's/^Chip type:[[:space:]]*//p' "${flash_log}" | head -1)"
}

# The NUCLEO-H753ZI's power-on prompt watches the ST-LINK; here the readers are
# already attached to the telemetry adapter, which stays enumerated, so the mark
# is taken when the board's own bridge goes and the report keeps its head.
await_esp32_power_cycle()
{
    local node="$1" waited=0
    echo ">>> UNPLUG the ${board_label}'s own USB cable now, leave it out for a few seconds, then plug it back in." >&2
    while [ "${waited}" -lt $((plug_timeout * 10)) ] && [ -e "${node}" ]; do
        sleep 0.1
        waited=$((waited + 1))
    done
    [ ! -e "${node}" ] \
        || fail 3 "power cycle timed out: ${node} did not disappear within ${plug_timeout} s, so the cable was not unplugged."
    stream_mark="$(wc -c < "${raw_stream}")"
    console_mark="$(wc -c < "${console_raw}")"
    waited=0
    while [ "${waited}" -lt $((plug_timeout * 10)) ] && [ ! -e "${node}" ]; do
        sleep 0.1
        waited=$((waited + 1))
    done
    [ -e "${node}" ] \
        || fail 3 "power cycle timed out: ${node} did not reappear within ${plug_timeout} s after the unplug."
    echo ">>> device node back: ${node}" >&2
}

# While the board is unpowered its TX line idles low, which the telemetry adapter
# reads as NUL bytes; they can land at the head of the build id line. The report
# is text and never carries a NUL, so they are dropped.
stream_since_mark()
{
    tail -c +"$((stream_mark + 1))" "${raw_stream}" 2>/dev/null | tr -d '\000'
}

# A closing line counts only after a build id line, since a reader attached to a
# running board can join a report midway. awk reads to the end rather than
# quitting at the first match, so tail is not killed by a closed pipe, which
# pipefail would report as a failure.
report_closed()
{
    stream_since_mark | awk -v marker="${report_end_marker}" '/^\[meta\] build_id=/ { opened = 1 } opened && index($0, marker) { seen = 1 } END { exit !seen }'
}

wait_for_report_esp32()
{
    local waited=0
    while [ "${waited}" -lt "${read_timeout}" ]; do
        [ "$(stream_since_mark | wc -c)" -lt "${byte_cap}" ] || return 0
        if report_closed; then
            return 0
        fi
        sleep 1
        waited=$((waited + 1))
    done
    return 0
}

# The first report after the mark, from its build id line through its closing
# verdict: the board reports on every boot, so the next one is cut away, and
# what came before the build id line is residue from before the restart.
cut_report()
{
    local all="${transcript}.all" start=""
    stream_since_mark | tr -d '\r' > "${all}"
    start="$(grep -n -m 1 '^\[meta\] build_id=' "${all}" | cut -d: -f1 || true)"
    discarded=0
    if [ -n "${start}" ]; then
        discarded="$(head -n "$((start - 1))" "${all}" | wc -c)"
    fi
    awk -v start="${start:-1}" -v marker="${report_end_marker}" 'NR >= start { print } NR >= start && index($0, marker) { exit }' "${all}" > "${transcript}"
    rm -f -- "${all}"
}

# Written before the refusals run, so a refused capture still leaves the
# console, where a panic or a probe overflow is reported.
write_console_record()
{
    {
        echo "# ctrlpp ${board_label} console, recorded beside $(basename "${destination}")"
        echo "# port:     $1"
        echo "# The board's own console: the framework's boot log, the stack probe and the"
        echo "# panic it ends in. Transported, not judged; the capture's refusals read only"
        echo "# the telemetry report."
        echo
        tail -c +"$((console_mark + 1))" "${console_raw}" | tr -d '\r'
    } > "${console_file}"
}

do_capture_esp32()
{
    [ -n "${destination}" ] || fail 2 "capture needs a destination file."
    [ -d "$(dirname "${destination}")" ] \
        || fail 2 "the destination directory '$(dirname "${destination}")' does not exist."
    [ -f "${bin_path}" ] || fail 3 "no ${bin_path}; run 'board_capture.sh --board=esp32 prepare' first."
    command -v python3 > /dev/null || fail 3 "python3 is needed to drive the board's reset lines."
    load_elf_reference

    local telemetry="" link="" console="" seconds=0
    telemetry="$(resolve_port)" || fail 3 "no telemetry port found: no /dev/serial/by-id/*FTDI* link and no CAPTURE_PORT."
    [ -c "${telemetry}" ] || fail 3 "'${telemetry}' is not a character device."
    link="$(console_link)" || fail 3 "no console port found: no /dev/serial/by-id/*1a86* or *CP210* link and no CAPTURE_CONSOLE_PORT."
    console="$(readlink -f "${link}")"
    [ -c "${console}" ] || fail 3 "'${console}' is not a character device."
    console_file="${destination%.txt}-console.txt"

    raw_stream="$(mktemp "${TMPDIR:-/tmp}/ctrlpp-board-capture-raw.XXXXXX")"
    transcript="$(mktemp "${TMPDIR:-/tmp}/ctrlpp-board-capture.XXXXXX")"
    console_raw="$(mktemp "${TMPDIR:-/tmp}/ctrlpp-board-console.XXXXXX")"
    flash_log="$(mktemp "${TMPDIR:-/tmp}/ctrlpp-board-flash.XXXXXX")"
    trap cleanup EXIT

    if [ "${capture_mode}" = "flash" ]; then
        check_descriptor "${bin_path}" "${reference_id}"
        flash_esp32 "${console}"
    fi

    seconds=$((read_timeout + console_tail + 5))
    [ "${capture_mode}" != "power-cycle" ] || seconds=$((seconds + 2 * plug_timeout))
    esp32_reader "${console}" "${console_raw}" "${seconds}"
    console_pid=$!
    esp32_reader "${telemetry}" "${raw_stream}" "${seconds}"
    reader_pid=$!

    case "${capture_mode}" in
        flash|reset-only)
            drive_en "${console}" hold
            sleep 1
            stream_mark="$(wc -c < "${raw_stream}")"
            console_mark="$(wc -c < "${console_raw}")"
            drive_en "${console}" release
            ;;
        power-cycle)
            await_esp32_power_cycle "${link}"
            console_stop
            console="$(readlink -f "${link}")"
            esp32_reader "${console}" "${console_raw}" "${seconds}"
            console_pid=$!
            ;;
    esac

    wait_for_report_esp32
    reader_stop
    if report_closed; then
        sleep "${console_tail}"
    fi
    console_stop

    write_console_record "${console}"
    cut_report
    verify_transcript
    write_artifact_esp32 "${telemetry}" "${console}"

    echo "capture written: ${destination}"
    echo "console written: ${console_file}"
    echo "build id:        ${reference_id}"
}

reset_note()
{
    case "${capture_mode}" in
        flash|reset-only)
            echo "# This run was started by this script's EN-pin pulse through the bridge's RTS"
            echo "# line. It was NOT a power cycle: the part reports an EN-pin reset and a power"
            echo "# cycle as the same power_on class, so the label cannot say which happened,"
            echo "# and neither this run nor a series of such runs is power-on evidence."
            ;;
        power-cycle)
            echo "# The operator was asked to unplug and replug the board, and its console node"
            echo "# vanished and returned. That proves the USB bridge lost its connection, not"
            echo "# that the supply decayed, and the label cannot tell the two apart either."
            ;;
        *)
            echo "# This script did not restart the board; the report is the first to begin"
            echo "# after the readers attached, and the label is the image's own."
            ;;
    esac
}

schedule_field()
{
    grep -m 1 '^\[meta\] schedule ' "${transcript}" | sed -n "s/.* $1=\([^ ]*\).*/\1/p"
}

# Assembled from the schedule and interruption lines, so it names only causes
# the posture reports and quotes the tick count rather than asserting what varied.
band_cause_esp32()
{
    local crossed=""
    crossed="$(sed -n 's/^\[verdict\] interruption region=\([^ ]*\) regime=warm tick_crossed=\([0-9]*\) of=\([0-9]*\)$/\1 \2 of \3/p' "${transcript}" | paste -sd ';' - | sed 's/;/; /g')"
    echo "#"
    echo "# Band cause: the timed code runs from $(schedule_field code), fetched through the"
    echo "# cache at flash_mode=$(schedule_field flash_mode) flash_freq=$(schedule_field flash_freq), so the cold-to-warm"
    echo "# difference includes cache fill. Within the warm series the spread includes"
    echo "# interrupts landing in a sample, with the scheduler ticking at $(schedule_field tick_hz) Hz and"
    echo "# the enabled-interrupt mask $(schedule_field intr_enabled) on the timed core, and the other core"
    echo "# competing for the same flash interface; this capture does not separate them."
    echo "# The tick landed in, per warm region: ${crossed:-no interruption line}."
}

write_artifact_esp32()
{
    {
        echo "# ctrlpp ${board_label} on-target capture"
        echo "# captured:  $(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo "# mode:      ${capture_mode}"
        echo "# telemetry: $1"
        echo "# console:   $2 (recorded in $(basename "${console_file}"))"
        echo "# image:     ${elf_path}"
        echo "# kind:      ${image_kind}"
        echo "# build id:  ${reference_id}"
        if [ "${image_kind}" = "timing" ]; then
            echo "# evidence:  ${reference_evidence_id} (the allocation-evidence image built from the same sources)"
        fi
        echo "# toolchain: ${reference_toolchain}"
        echo "# framework: $(framework_line), as the build recorded it"
        echo "# chip:      ${chip_line:-not read in this mode; esptool reports it when flashing}"
        echo "#"
        echo "# The build id above is the sha256 digest of the ELF, which the framework's"
        echo "# build writes into the image descriptor and the image prints; sha256sum of"
        echo "# the ELF reproduces it. It is an identity check against accidental"
        echo "# staleness, not a signature."
        if [ "${capture_mode}" = "flash" ]; then
            echo "# Before flashing, the descriptor of $(basename "${bin_path}") was read with esptool"
            echo "# and named the same digest, and esptool verified every region it wrote."
        fi
        echo "#"
        echo "# Reset: $(posture_of "${transcript}" | grep -o "${reset_quote}" || true)."
        reset_note
        echo "#"
        echo "# The [exact] rows print each float as its binary32 bit pattern in 8"
        echo "# uppercase hexadecimal digits, so two runs that differ in the last bit differ"
        echo "# in the row."
        echo "#"
        echo "# Stack: the [verdict] stack line is the control task's high-water as FreeRTOS"
        echo "# counts it, against the ${reference_reserve}-byte task stack the source names, taken"
        echo "# with $(posture_of "${transcript}" | grep -o 'stack_watchpoint=[^ ]*' || true). The evidence image's stack probe prints its"
        echo "# per-dimension lines in the console record, and they are not judged here."
        if [ "${image_kind}" = "timing" ]; then
            band_cause_esp32
        fi
        echo "#"
        echo "# Report window: ${discarded} bytes before the report's build id line were"
        echo "# discarded, the previous boot's residue or line noise across the restart."
        echo
        cat "${transcript}"
    } > "${destination}"
}

# --- capture ------------------------------------------------------------------

do_capture()
{
    if [ "${board}" = "esp32" ]; then
        do_capture_esp32
        return 0
    fi
    [ -n "${destination}" ] || fail 2 "capture needs a destination file."
    [ -d "$(dirname "${destination}")" ] \
        || fail 2 "the destination directory '$(dirname "${destination}")' does not exist."
    [ -f "${bin_path}" ] || fail 3 "no ${bin_path}; run 'board_capture.sh prepare' first."
    load_elf_reference

    local port

    port="$(resolve_port)" || fail 3 "no ST-LINK serial port found: neither /dev/serial/by-id/*STLINK* nor a USB device 0483:374e with a tty child resolved."
    [ -c "${port}" ] || fail 3 "'${port}' is not a character device."

    if [ "${capture_mode}" = "flash" ]; then
        [ -d "${mount_dir}" ] || fail 3 "'${mount_dir}' is not a directory; the ST-LINK mass storage is not mounted."
        mountpoint -q "${mount_dir}" || fail 3 "'${mount_dir}' exists but is not a mount point; refusing to write the image into the local filesystem."
    fi

    load_handshake_byte
    raw_stream="$(mktemp "${TMPDIR:-/tmp}/ctrlpp-board-capture-raw.XXXXXX")"
    transcript="$(mktemp "${TMPDIR:-/tmp}/ctrlpp-board-capture.XXXXXX")"
    trap cleanup EXIT

    reader_start "${port}"

    case "${capture_mode}" in
        flash)
            cp -- "${bin_path}" "${mount_dir}/" || fail 3 "the copy onto ${mount_dir} failed."
            sync
            wait_for_port_cycle "${port}"
            reader_stop
            reader_start "${port}"
            ;;
        reset-only)
            board_reset
            ;;
        no-reset)
            ;;
    esac

    wait_for_report "${port}"
    reader_stop

    tr -d '\r' < "${raw_stream}" > "${transcript}"
    verify_transcript
    write_artifact "${port}" "${reference_id}"

    echo "capture written: ${destination}"
    echo "build id:        ${reference_id}"
}

# The same refusals over a file already on disk. A capture file's header lines
# start with '#', which no refusal's anchored pattern matches.
do_check()
{
    [ -n "${destination}" ] || fail 2 "check needs a transcript file."
    [ -f "${destination}" ] || fail 2 "'${destination}' is not a file."
    load_elf_reference

    check_file "${destination}"
    [ "${check_status}" -eq 0 ] || exit "${check_status}"
    echo "check passed: ${destination}"
}

# Every refusal of the chosen mode over one file, against the loaded reference.
# A refusal exits, so the verification runs in a subshell and its status lands
# in check_status. The function itself always returns 0: called on the left of
# || it would suspend errexit inside the subshell too, and a verification that
# runs on past a failed command is not one.
check_status=0

check_file()
{
    transcript="$(mktemp "${TMPDIR:-/tmp}/ctrlpp-board-check.XXXXXX")"
    tr -d '\r' < "$1" > "${transcript}"
    set +e
    ( set -e; verify_transcript )
    check_status=$?
    set -e
    rm -f -- "${transcript}"
    transcript=""
    return 0
}

# --- compare ------------------------------------------------------------------
#
# A series is N restarts of one resident image, each captured and refused on its
# own, then compared: the exact-class block byte for byte, the spread-class block
# field by field against a stated bound. The capture subcommand does the
# capturing, so every run of a series passes every refusal a single capture does.

script_path="${script_dir}/$(basename "${BASH_SOURCE[0]}")"
compare_attempts="${COMPARE_ATTEMPTS:-3}"
max_refused="${COMPARE_MAX_REFUSED:-}"
plug_timeout="${COMPARE_PLUG_TIMEOUT:-120}"

# Every spread-class column belongs to one of these. A shape field says what was
# measured, so runs that disagree on it measured different things and are not
# compared; a bounded field is a cycle count held to the bound; a reported field
# is published with its spread but is not a cycle count the bound can apply to.
spread_key_fields="region regime"
spread_shape_fields="n n_basis block wraps"
spread_bounded_fields="raw_min raw_median raw_max"
spread_reported_fields="max_at per_iteration_raw_max"

requested_class=""
decoded_class=""
series_count=""
series_dir=""
compare_tmp=""
compare_status=0
compare_note=""
runs=()
run_log=()

# The first failure decides the exit status; the comparison still completes so
# the summary reports every spread, including on a series that failed.
compare_fail()
{
    echo "REFUSED: $2" >&2
    if [ "${compare_status}" -eq 0 ]; then
        compare_status="$1"
        compare_note="$2"
    fi
    return 0
}

# Anchored on the tag, so the '# Reset:' header line of a capture file is not
# mistaken for the image's own report.
posture_of()
{
    tr -d '\r' < "$1" | grep -m 1 '^\[meta\] posture ' || true
}

posture_field()
{
    posture_of "$1" | sed -n "s/.* $2=\([^ ]*\).*/\1/p"
}

block_of()
{
    tr -d '\r' < "$2" | grep "^\[$1\] " || true
}

# The node the operator's replug makes vanish and return. The by-id link is
# stable across a re-enumeration where the tty number is not.
port_node()
{
    local link
    if [ -n "${CAPTURE_PORT:-}" ]; then
        printf '%s\n' "${CAPTURE_PORT}"
        return 0
    fi
    for link in /dev/serial/by-id/*STLINK*; do
        if [ -e "${link}" ]; then
            printf '%s\n' "${link}"
            return 0
        fi
    done
    resolve_port
}

# The script cannot cut the target's power, so it asks for the cycle and waits,
# bounded both ways, for the node to go and come back. The prompt is only an
# instruction: the decoded class on the run's posture line is what counts.
await_power_cycle()
{
    local node="" waited=0
    node="$(port_node)" || fail 3 "no ST-LINK serial port found to watch for the power cycle."
    echo ">>> UNPLUG the board's ST-LINK USB cable now, leave it out for a few seconds, then plug it back in." >&2
    while [ "${waited}" -lt $((plug_timeout * 10)) ] && [ -e "${node}" ]; do
        sleep 0.1
        waited=$((waited + 1))
    done
    [ ! -e "${node}" ] \
        || fail 3 "power cycle timed out: ${node} did not disappear within ${plug_timeout} s, so the cable was not unplugged."
    waited=0
    while [ "${waited}" -lt $((plug_timeout * 10)) ] && [ ! -e "${node}" ]; do
        sleep 0.1
        waited=$((waited + 1))
    done
    [ -e "${node}" ] \
        || fail 3 "power cycle timed out: ${node} did not reappear within ${plug_timeout} s after the unplug."
    echo ">>> device node back: ${node}" >&2
    await_node_open "${node}"
}

# A node can appear before the device manager has made it accessible, so it is
# opened once, bounded, before the capture's reader needs it. The image waits
# for the handshake meanwhile, so the wait costs no part of the report.
await_node_open()
{
    local waited=0
    until stty -F "$1" 115200 raw -echo 2>/dev/null; do
        [ "${waited}" -lt $((plug_timeout * 10)) ] \
            || fail 3 "power cycle: $1 came back but could not be opened within ${plug_timeout} s."
        sleep 0.1
        waited=$((waited + 1))
    done
}

# A lost report head (an empty stream, or no fingerprint line) is the transport
# failing, not the run, so it is retried; any other refusal ends the series.
collect_one()
{
    local file="$1" mode="$2" attempt=1 status=0
    while [ "${attempt}" -le "${compare_attempts}" ]; do
        status=0
        if [ "${mode}" = "power_on" ] && [ "${board}" = "esp32" ]; then
            "${script_path}" capture ${board_flag} ${kind_flag} --power-cycle "${file}" || status=$?
        elif [ "${mode}" = "power_on" ]; then
            await_power_cycle
            "${script_path}" capture ${kind_flag} --no-reset "${file}" || status=$?
        else
            "${script_path}" capture ${board_flag} ${kind_flag} --reset-only "${file}" || status=$?
        fi
        if [ "${status}" -ne 4 ] && [ "${status}" -ne 5 ]; then
            run_log+=("$(basename "${file}") attempts=${attempt} capture_status=${status}")
            return "${status}"
        fi
        echo "capture of $(basename "${file}") lost the report head (status ${status}); attempt ${attempt} of ${compare_attempts}." >&2
        attempt=$((attempt + 1))
    done
    run_log+=("$(basename "${file}") attempts=${compare_attempts} capture_status=${status}")
    return "${status}"
}

# The priming run is not counted. The image clears the reset flags after it
# reports them, so a run after it labels only its own reset; without it the first
# counted run could carry flags accumulated under an image that never cleared.
# The ESP32 rewrites its reset code on every reset and needs no clearing, but it
# is primed the same way, an EN-pin restart, so the resident image is checked
# before the operator is asked for anything.
#
# A run that decoded another class is evidence about the restart, not about the
# image, so it is kept under its own name, listed in the summary and collected
# again; a run the capture refused ends the series, and the summary is still
# written for the runs that were counted.
collect_series()
{
    local index=1 refused=0 file="" status=0 class="" kept=""
    collect_one "${series_dir}/priming.txt" pin || fail "$?" "the priming run was refused; the series was not started."
    while [ "${index}" -le "${series_count}" ]; do
        file="$(printf '%s/run-%02d.txt' "${series_dir}" "${index}")"
        status=0
        collect_one "${file}" "${requested_class}" || status=$?
        if [ "${status}" -ne 0 ]; then
            compare_fail "${status}" "run ${index} of the series was refused by the capture; the series stops at $((index - 1)) counted runs."
            return 0
        fi
        class="$(posture_field "${file}" reset)"
        if [ "${class}" = "${decoded_class}" ]; then
            index=$((index + 1))
            continue
        fi
        refused=$((refused + 1))
        kept="$(printf '%s/refused-%02d.txt' "${series_dir}" "${refused}")"
        mv -- "${file}" "${kept}"
        run_log[${#run_log[@]} - 1]="$(basename "${kept}") (collected as $(basename "${file}")) ${run_log[${#run_log[@]} - 1]#* }"
        if [ -n "${max_refused}" ] && [ "${refused}" -ge "${max_refused}" ]; then
            compare_fail 26 "run ${index} decoded reset=${class:-absent}, and ${refused} runs of another class reach COMPARE_MAX_REFUSED; the series stops at $((index - 1)) counted runs."
            return 0
        fi
        echo "run ${index} decoded reset=${class:-absent} where the series requires reset=${decoded_class}; it is kept as $(basename "${kept}"), not counted, and run ${index} is collected again." >&2
    done
    return 0
}

check_classes()
{
    local run="" class=""
    for run in "${runs[@]}"; do
        class="$(posture_field "${run}" reset)"
        [ "${class}" = "${decoded_class}" ] \
            || compare_fail 26 "$(basename "${run}") decoded reset=${class:-absent} where the series requires reset=${decoded_class}."
    done
    return 0
}

# The ESP32 decodes an EN-pin pulse and a power cycle alike, so the label cannot
# keep the two series apart; the capture mode each run's header records does.
restart_mode()
{
    case "$1" in
        en_pin) echo "reset-only" ;;
        power_on) echo "power-cycle" ;;
    esac
}

check_restarts()
{
    local run="" mode="" want=""
    want="$(restart_mode "${requested_class}")"
    [ "${board}" = "esp32" ] && [ -n "${want}" ] || return 0
    for run in "${runs[@]}"; do
        mode="$(sed -n 's/^# mode:[[:space:]]*//p' "${run}" | head -1)"
        [ "${mode}" = "${want}" ] \
            || compare_fail 26 "$(basename "${run}") records the restart mode '${mode:-none}' where this ${requested_class} series requires '${want}'; both decode reset=power_on, so the runs of one are not pooled with the other."
    done
    return 0
}

restart_line()
{
    case "${requested_class}" in
        en_pin)
            echo "# restart:       this script's EN-pin pulse through the console bridge's RTS line, NOT a power"
            echo "#                cycle. The image decodes both as reset=power_on, so the decoded label cannot"
            echo "#                distinguish them, and this series is not power-on evidence."
            ;;
        power_on)
            echo "# restart:       a power cycle the operator was asked for; the decoded class is the only"
            echo "#                evidence of it$([ "${board}" = "esp32" ] && echo ", and on this leg it cannot tell a power cycle from an EN-pin reset")."
            ;;
        pin)
            echo "# restart:       this script's st-flash reset."
            ;;
        *)
            echo "# restart:       not recorded by this comparison; the decoded class is the image's own."
            ;;
    esac
}

# The fingerprint and every posture field but the reset's own must agree, since
# a spread is only reported against the one posture that produced it.
check_shared()
{
    local run="" first_id="" id="" first_posture="" posture=""
    first_id="$(block_of meta "${runs[0]}" | sed -n 's/^\[meta\] build_id=\([0-9a-f]*\).*/\1/p' | head -1)"
    first_posture="$(posture_of "${runs[0]}" | sed -E 's/ reset(_flags|_rsr|_reason|_rom)?=[^ ]*//g')"
    for run in "${runs[@]}"; do
        id="$(block_of meta "${run}" | sed -n 's/^\[meta\] build_id=\([0-9a-f]*\).*/\1/p' | head -1)"
        [ "${id}" = "${first_id}" ] \
            || compare_fail 29 "$(basename "${run}") carries build id ${id:-absent} where $(basename "${runs[0]}") carries ${first_id:-absent}; a series is one image."
        posture="$(posture_of "${run}" | sed -E 's/ reset(_flags|_rsr|_reason|_rom)?=[^ ]*//g')"
        [ "${posture}" = "${first_posture}" ] \
            || compare_fail 29 "$(basename "${run}") reports the posture '${posture}' where $(basename "${runs[0]}") reports '${first_posture}'."
    done
    shared_id="${first_id}"
    shared_posture="${first_posture}"
    return 0
}

# Byte for byte and positional: no rounding, no normalization, no sorting, so a
# one-bit difference or a reordered row is a difference.
compare_exact()
{
    local run="" where="" line="" reference="${compare_tmp}/exact-reference" candidate="${compare_tmp}/exact-candidate"
    block_of exact "${runs[0]}" > "${reference}"
    exact_rows=$(($(wc -l < "${reference}") - 1))
    exact_digest="$(sha256sum < "${reference}" | cut -d' ' -f1)"
    exact_verdict="byte-identical"
    for run in "${runs[@]:1}"; do
        block_of exact "${run}" > "${candidate}"
        where="$(cmp -- "${reference}" "${candidate}" 2>&1)" && continue
        line="$(printf '%s\n' "${where}" | sed -n 's/.*line \([0-9][0-9]*\).*/\1/p')"
        exact_verdict="DIFFERS"
        if printf '%s\n' "${where}" | grep -q 'EOF on'; then
            compare_fail 27 "$(basename "${run}") and $(basename "${runs[0]}") carry exact-class blocks of different lengths: ${where##*: }."
            continue
        fi
        compare_fail 27 "$(basename "${run}") differs from $(basename "${runs[0]}") in the exact-class block at block line ${line:-?}, byte $(first_differing_byte "${reference}" "${candidate}" "${line:-1}") of that line: '$(sed -n "${line:-1}p" "${reference}")' against '$(sed -n "${line:-1}p" "${candidate}")'."
    done
    return 0
}

first_differing_byte()
{
    awk -v n="$3" 'FNR == n { line[++seen] = $0 } END {
        a = line[1]; b = line[2]; i = 1
        while (i <= length(a) && substr(a, i, 1) == substr(b, i, 1)) ++i
        print i }' "$1" "$2"
}

# One table row per spread row and field: the minimum, maximum and spread across
# the runs, the bound where the field is a cycle count, and the verdict.
compare_spread()
{
    local index=0 run="" files=() problem=""
    for run in "${runs[@]}"; do
        index=$((index + 1))
        block_of spread "${run}" | sed 's/^\[spread\] //' > "${compare_tmp}/spread-${index}"
        files+=("${compare_tmp}/spread-${index}")
    done
    awk -F, -v key="${spread_key_fields}" -v shape="${spread_shape_fields}" \
        -v bounded="${spread_bounded_fields}" -v reported="${spread_reported_fields}" \
        -v bound="${spread_bound:-}" -f /dev/stdin "${files[@]}" > "${compare_tmp}/spread-table" <<'AWK'
function classify(list, kind,    n, i, names) { n = split(list, names, " "); for (i = 1; i <= n; ++i) cls[names[i]] = kind }
function problem(status, text) { if (!reported_problem) print "PROBLEM\t" status "\t" text; reported_problem = 1 }
function decimals(v) { return index(v, ".") ? length(v) - index(v, ".") : 0 }
BEGIN { classify(key, "key"); classify(shape, "shape"); classify(bounded, "bounded"); classify(reported, "reported") }
FNR == 1 {
    if (run > 0 && count[run] != count[1]) problem(29, "run " run " carries " count[run] " spread rows where run 1 carries " count[1] ".")
    ++run
    if (run == 1) {
        header = $0; width = NF
        for (i = 1; i <= NF; ++i) { name[i] = $i; if (!($i in cls)) problem(29, "the spread-class column '" $i "' belongs to no class this comparison knows, so it can be neither bounded nor required equal.") }
    } else if ($0 != header) problem(29, "run " run " carries the spread header '" $0 "' where run 1 carries '" header "'.")
    next
}
{
    row = ++count[run]
    if (run == 1) keys[row] = $1 "," $2
    else if (keys[row] != $1 "," $2) problem(29, "run " run " carries the spread row '" $1 "," $2 "' at position " row " where run 1 carries '" keys[row] "'.")
    for (i = 3; i <= width; ++i) {
        kind = cls[name[i]]
        if (kind == "shape") {
            if (run == 1) first[row, i] = $i
            else if ($i != first[row, i]) problem(29, "run " run " measured " keys[row] " with " name[i] "=" $i " where run 1 used " first[row, i] ".")
            continue
        }
        if ($i !~ /^[0-9]+(\.[0-9]+)?$/) { problem(25, "run " run " carries the non-numeric " name[i] " '" $i "' in the spread row " keys[row] "."); continue }
        if (run == 1 || $i + 0 < lo[row, i] + 0) lo[row, i] = $i
        if (run == 1 || $i + 0 > hi[row, i] + 0) hi[row, i] = $i
    }
}
END {
    if (run > 1 && count[run] != count[1]) problem(29, "run " run " carries " count[run] " spread rows where run 1 carries " count[1] ".")
    for (row = 1; row <= count[1]; ++row)
        for (i = 3; i <= width; ++i) {
            kind = cls[name[i]]
            if (kind != "bounded" && kind != "reported") continue
            d = decimals(lo[row, i]) > decimals(hi[row, i]) ? decimals(lo[row, i]) : decimals(hi[row, i])
            spread = sprintf("%." d "f", hi[row, i] - lo[row, i])
            verdict = "reported"
            limit = "-"
            if (kind == "bounded") {
                limit = bound
                verdict = (spread + 0 <= bound + 0) ? "within" : "EXCEEDS"
                if (verdict == "EXCEEDS") problem(28, "the spread of " name[i] " for " keys[row] " is " spread " cycles across " run " runs, past the stated bound of " bound " cycles (" lo[row, i] " to " hi[row, i] ").")
            }
            print "ROW\t" keys[row] "," name[i] "," lo[row, i] "," hi[row, i] "," spread "," limit "," verdict
        }
}
AWK
    problem="$(grep -m 1 '^PROBLEM' "${compare_tmp}/spread-table" || true)"
    if [ -n "${problem}" ]; then
        compare_fail "$(printf '%s\n' "${problem}" | cut -f2)" "$(printf '%s\n' "${problem}" | cut -f3-)"
    fi
    spread_rows="$(sed -n 's/^ROW\t//p' "${compare_tmp}/spread-table")"
    return 0
}

write_summary()
{
    local run="" verdict="PASS" refused=()
    for run in "${series_dir}"/refused-[0-9][0-9].txt; do
        [ -f "${run}" ] && refused+=("${run}")
    done
    [ "${compare_status}" -eq 0 ] || verdict="FAIL (status ${compare_status}): ${compare_note}"
    {
        echo "# ctrlpp ${board_label} determinism comparison"
        echo "# compared:      $(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo "# mode:          $([ "${replay}" = "yes" ] && echo "replay of transcripts on disk" || echo "collected by this invocation")"
        echo "# kind:          ${image_kind}"
        echo "# reset class:   ${requested_class} (requested; each run must decode reset=${decoded_class}; every counted and every refused run's decoded class is listed below)"
        restart_line
        echo "# runs:          ${#runs[@]} counted$([ -n "${series_count}" ] && echo " of ${series_count} requested")"
        echo "# refused:       ${#refused[@]} decoded another class and were not counted"
        echo "# build id:      ${shared_id:-absent} (shared by every counted run)"
        if [ "${image_kind}" = "timing" ]; then
            echo "# evidence:      ${reference_evidence_id:-absent} (the allocation-evidence image built from the same sources)"
        fi
        echo "# toolchain:     ${reference_toolchain:-absent}"
        echo "# posture:       ${shared_posture:-absent}"
        echo "# host load:     $(cut -d' ' -f1-3 /proc/loadavg 2>/dev/null || echo unknown) (the figures count core cycles on the target, not host time)"
        echo "#"
        echo "# Per run, as each image decoded its own reset status register:"
        if [ -f "${series_dir}/priming.txt" ] && [ "${replay}" = "no" ]; then
            echo "#   priming.txt (not counted) $(posture_of "${series_dir}/priming.txt" | grep -o "${reset_quote}")"
        fi
        for run in "${runs[@]}"; do
            echo "#   $(basename "${run}") $(posture_of "${run}" | grep -o "${reset_quote}")"
        done
        for run in "${refused[@]+"${refused[@]}"}"; do
            echo "#   $(basename "${run}") (refused, not counted) $(posture_of "${run}" | grep -o "${reset_quote}" || echo "no posture line")"
        done
        for run in "${run_log[@]+"${run_log[@]}"}"; do
            echo "#   collection: ${run}"
        done
        echo "#"
        echo "# Exact class: ${exact_rows:-0} rows plus the header, ${exact_verdict:-not compared} across ${#runs[@]} runs."
        echo "#   sha256 of the block in $(basename "${runs[0]:-none}"): ${exact_digest:-none}"
        echo "#"
        echo "# Spread class: bound ${spread_bound:-none} cycles on ${spread_bounded_fields}; ${spread_reported_fields} are"
        echo "# reported without a bound; ${spread_shape_fields} must agree across runs."
        echo "region,regime,field,min,max,spread,bound,verdict"
        if [ -n "${spread_rows:-}" ]; then
            printf '%s\n' "${spread_rows}"
        else
            echo "# the spread-class block carries no rows on this image"
        fi
        echo "#"
        echo "# Verdict: ${verdict}"
    } > "${series_dir}/compare-summary.txt"
}

refuse_en_pin_off_esp32()
{
    [ "${requested_class}" != "en_pin" ] || [ "${board}" = "esp32" ] \
        || { usage; fail 2 "en_pin is the ESP32 leg's EN-pin restart class; the ${board_label} has no EN pin this script drives, and its script-driven restart is the 'pin' class."; }
}

compare_arguments()
{
    if [ "${replay}" = "yes" ]; then
        [ "${#operands[@]}" -eq 2 ] || { usage; fail 2 "compare --replay takes a reset class and a series directory."; }
        requested_class="${operands[0]}"
        series_dir="${operands[1]}"
        refuse_en_pin_off_esp32
    else
        [ "${#operands[@]}" -eq 3 ] || { usage; fail 2 "compare takes a run count, a reset class and a destination directory."; }
        series_count="${operands[0]}"
        requested_class="${operands[1]}"
        series_dir="${operands[2]}"
        refuse_en_pin_off_esp32
        printf '%s\n' "${series_count}" | grep -Eqx '[1-9][0-9]*' && [ "${series_count}" -ge 2 ] \
            || fail 2 "the run count must be an integer of at least 2, since one run compares with nothing: '${series_count}'."
        case " ${collectable_classes} " in
            *" ${requested_class} "*) ;;
            *) fail 2 "on the ${board_label} this station can collect a series of these classes only: ${collectable_classes}; not '${requested_class}'. A series of another class can be compared with --replay." ;;
        esac
        printf '%s\n' "${plug_timeout}" | grep -Eqx '[1-9][0-9]*' \
            || fail 2 "COMPARE_PLUG_TIMEOUT must be a whole number of seconds: '${plug_timeout}'."
        printf '%s\n' "${compare_attempts}" | grep -Eqx '[1-9][0-9]*' \
            || fail 2 "COMPARE_ATTEMPTS must be a positive integer: '${compare_attempts}'."
        [ -z "${max_refused}" ] || printf '%s\n' "${max_refused}" | grep -Eqx '[1-9][0-9]*' \
            || fail 2 "COMPARE_MAX_REFUSED must be a positive integer: '${max_refused}'."
    fi
    printf '%s\n' "${requested_class}" | grep -Eqx '[a-z_]+' || fail 2 "'${requested_class}' is not a reset class name."
    decoded_class="${requested_class}"
    [ "${requested_class}" != "en_pin" ] || decoded_class="power_on"
    if [ "${image_kind}" = "timing" ]; then
        printf '%s\n' "${spread_bound}" | grep -Eqx '[0-9]+' \
            || fail 2 "a timing series needs its spread bound stated as --bound=<cycles>, a non-negative integer: '${spread_bound}'."
    else
        [ -z "${spread_bound}" ] || fail 2 "--bound applies to the timing image only; the evidence image times nothing."
    fi
    kind_flag=""
    [ "${image_kind}" = "timing" ] && kind_flag="--timing"
    board_flag=""
    [ "${board}" = "esp32" ] && board_flag="--board=esp32"
    return 0
}

do_compare()
{
    compare_arguments
    load_elf_reference
    run_compare
}

run_compare()
{
    compare_tmp="$(mktemp -d "${TMPDIR:-/tmp}/ctrlpp-board-compare.XXXXXX")"
    trap 'cleanup; rm -rf -- "${compare_tmp}"' EXIT

    local run=""
    if [ "${replay}" = "yes" ]; then
        [ -d "${series_dir}" ] || fail 2 "'${series_dir}' is not a directory."
        for run in "${series_dir}"/run-[0-9][0-9].txt; do
            [ -f "${run}" ] || continue
            check_file "${run}"
            [ "${check_status}" -eq 0 ] || fail "${check_status}" "$(basename "${run}") fails a single-capture refusal; it is not compared."
            runs+=("${run}")
        done
        [ "${#runs[@]}" -ge 2 ] || fail 2 "'${series_dir}' holds ${#runs[@]} run-NN.txt files; a comparison needs at least 2."
    else
        mkdir -p -- "${series_dir}" || fail 2 "cannot create '${series_dir}'."
        ! compgen -G "${series_dir}/run-*.txt" > /dev/null && ! compgen -G "${series_dir}/refused-*.txt" > /dev/null \
            || fail 2 "'${series_dir}' already holds run files; a series is collected into a directory of its own."
        collect_series
        for run in "${series_dir}"/run-[0-9][0-9].txt; do
            [ -f "${run}" ] && runs+=("${run}")
        done
    fi

    check_classes
    check_restarts
    if [ "${#runs[@]}" -gt 0 ]; then
        check_shared
        compare_exact
        compare_spread
    fi
    write_summary
    echo "summary written: ${series_dir}/compare-summary.txt"
    [ "${compare_status}" -eq 0 ] || exit "${compare_status}"
    if [ "${requested_class}" = "en_pin" ]; then
        echo "series PASS: ${#runs[@]} runs, each restarted by an EN-pin pulse (not a power cycle) and decoded reset=power_on, build id ${shared_id}"
    else
        echo "series PASS: ${#runs[@]} runs, reset=${requested_class}, build id ${shared_id}"
    fi
}

# --- discover -----------------------------------------------------------------

do_discover()
{
    local port="" console=""
    if [ "${board}" = "esp32" ]; then
        if port="$(resolve_port)" && [ -c "${port}" ] && console="$(resolve_console_port)" && [ -c "${console}" ]; then
            echo "board found: console ${console}, telemetry ${port}. Nothing was flashed or read here, so no verdict is formed; 'capture' takes the evidence."
            return 0
        fi
        echo "SKIP: no ${board_label} found (no console bridge under /dev/serial/by-id/*1a86* or *CP210* and telemetry adapter under *FTDI*, and no CAPTURE_CONSOLE_PORT and CAPTURE_PORT naming character devices); the on-target check was not performed: no image was flashed, no boot report was read and no resident image was fingerprinted."
        return 0
    fi
    if port="$(resolve_port)" && [ -c "${port}" ]; then
        echo "board found: ${port}. Nothing was flashed or read here, so no verdict is formed; 'capture' takes the evidence."
        return 0
    fi
    echo "SKIP: no board found (no /dev/serial/by-id/*STLINK* link, no USB device 0483:374e with a tty, no CAPTURE_PORT naming a character device); the on-target check was not performed: no image was flashed, no boot report was read and no resident image was fingerprinted."
}

# --- self-test ----------------------------------------------------------------
#
# Each red case is one corruption synthesized into a temporary copy of the
# committed recordings, asserted to draw one refusal status and its message.
# Corruptions are synthesized rather than committed, so covering a new refusal
# is one more line below and the fixtures stay the recordings themselves.

fixture_dir="${script_dir}/board_capture_fixtures"
self_test_root=""
self_test_cases=0
self_test_failures=0
self_test_green=0
self_test_red=()

fixture_value()
{
    sed -n "s/^$1=//p" "${self_test_root}/${2:-reference.txt}" | head -1
}

# A recording named esp32-<kind> is the ESP32 leg's and is checked against that
# leg's references. Every case runs in a subshell, so the board selected here
# does not outlive the case.
use_fixture_reference()
{
    local kind="$1" file="reference.txt"
    if [ "${kind#esp32-}" != "${kind}" ]; then
        select_board esp32
        kind="${kind#esp32-}"
        file="esp32-reference.txt"
    fi
    image_kind="${kind}"
    reference_id="$(fixture_value "${kind}_build_id" "${file}")"
    reference_reserve="$(fixture_value stack_reserve "${file}")"
    reference_evidence_id="$(fixture_value evidence_build_id "${file}")"
    reference_toolchain="$(fixture_value toolchain "${file}")"
    reference_predates="$(fixture_value predates "${file}")"
    reference_origin="board_capture_fixtures/${file}"
}

self_test_check()
{
    use_fixture_reference "$1"
    check_file "$2"
    return "${check_status}"
}

# The checks a built image gets: nothing is predated, so every line a built
# image prints is required of the recording.
self_test_check_current()
{
    use_fixture_reference "$1"
    reference_predates=""
    check_file "$2"
    return "${check_status}"
}

self_test_compare()
{
    use_fixture_reference "$1"
    replay="yes"
    operands=("$2" "$4")
    spread_bound="$3"
    compare_arguments
    run_compare
}

# A refusal exits, so each case runs in a subshell; errexit is restated there
# for the reason check_file gives.
expect_status()
{
    local want="$1" fragment="$2" label="$3" status=0 log="${self_test_root}/case.log" verdict="red"
    shift 3
    self_test_cases=$((self_test_cases + 1))
    set +e
    ( set -e; "$@" ) > /dev/null 2> "${log}"
    status=$?
    set -e
    if [ "${status}" -ne "${want}" ] || { [ -n "${fragment}" ] && ! grep -qF -- "${fragment}" "${log}"; }; then
        self_test_failures=$((self_test_failures + 1))
        printf '  FAILED %2d  %s: exited %d where %d with "%s" was expected\n' "${want}" "${label}" "${status}" "${want}" "${fragment}"
        sed 's/^/             /' "${log}"
        return 0
    fi
    if [ "${want}" -eq 0 ]; then
        verdict="green"
        self_test_green=$((self_test_green + 1))
    else
        self_test_red+=("${want}")
    fi
    printf '  %-6s %2d  %s\n' "${verdict}" "${want}" "${label}"
}

corrupt()
{
    local file="${self_test_root}/case-${self_test_cases}.txt"
    sed -e "$2" "${self_test_root}/$1.txt" > "${file}"
    printf '%s\n' "${file}"
}

# Three runs of one recording, the second corrupted. A third argument names the
# capture mode each run's header records, as a capture file's header would, so
# the corruption can reach the header too.
series()
{
    local dir="${self_test_root}/series-${self_test_cases}" index="" source="${self_test_root}/$1.txt"
    mkdir -- "${dir}"
    if [ -n "${3:-}" ]; then
        source="${dir}/source"
        { echo "# mode:      $3"; cat "${self_test_root}/$1.txt"; } > "${source}"
    fi
    for index in 01 03; do
        cp -- "${source}" "${dir}/run-${index}.txt"
    done
    sed -e "$2" "${source}" > "${dir}/run-02.txt"
    printf '%s\n' "${dir}"
}

red_check()
{
    local file=""
    file="$(corrupt "$1" "$5")"
    expect_status "$2" "$3" "$4" self_test_check "$1" "${file}"
}

red_current()
{
    local file=""
    file="$(corrupt "$1" "$5")"
    expect_status "$2" "$3" "$4" self_test_check_current "$1" "${file}"
}

red_series()
{
    local dir=""
    dir="$(series "$1" "$7" "${8:-}")"
    expect_status "$4" "$5" "$6" self_test_compare "$1" "$2" "$3" "${dir}"
}

# A replay of a series, then the summary line the case names, echoed to the log
# the case's fragment is matched against.
self_test_compare_summary()
{
    self_test_compare "$1" "$2" "$3" "$4"
    grep -F -- "$5" "$4/compare-summary.txt" >&2
}

# The ESP32 cut over a raw telemetry stream: an earlier report before the mark,
# then the NUL bytes an unpowered board's idle-low TX line reads as, ahead of
# the build id line, which the cut must still find.
self_test_cut()
{
    use_fixture_reference "$1"
    raw_stream="${self_test_root}/raw-${self_test_cases}"
    transcript="${self_test_root}/cut-${self_test_cases}"
    head -3 "${self_test_root}/$1.txt" > "${raw_stream}"
    stream_mark="$(wc -c < "${raw_stream}")"
    { printf '\000\000\000\000'; sed 's/$/\r/' "${self_test_root}/$1.txt"; cat "${self_test_root}/$1.txt"; } >> "${raw_stream}"
    report_closed
    cut_report
    [ "${discarded}" -eq 0 ]
    cp -- "${transcript}" "${transcript}.txt"
    check_file "${transcript}.txt"
    return "${check_status}"
}

# Read from this script's own source, so a refusal added without a case here
# leaves the self-test red rather than quietly uncovered.
implemented_refusals()
{
    {
        grep -oE '(^|[^[:alnum:]_])(compare_)?fail [0-9]+' "${script_path}" | grep -oE '[0-9]+$'
        grep -oE 'problem\([0-9]+,' "${script_path}" | grep -oE '[0-9]+'
    } | awk '$1 >= 4' | sort -n -u
}

self_test_cases_single()
{
    local reserve=""
    reserve="$(fixture_value stack_reserve)"
    red_check evidence 4 "EMPTY stream" "empty stream" 'd'
    red_check evidence 5 "carries no '[meta] build_id=' line" "fingerprint line missing" '/^\[meta\] build_id=/d'
    red_check evidence 6 "is NOT the built one" "fingerprint of another image" 's/^\(\[meta\] build_id=\)[0-9a-f]*/\10000000000000000000000000000000000000000/'
    red_check evidence 7 "verdict line -- it is truncated" "golden verdict missing" '/^\[verdict\] golden diff /d'
    red_check evidence 8 "golden verdict is FAIL" "golden verdict FAIL" 's/^\(\[verdict\] golden diff \)PASS/\1FAIL/'
    red_check evidence 9 "carries no '[canary] allocation' line" "allocation canary missing" '/^\[canary\] allocation /d'
    red_check evidence 10 "canary did not pass" "allocation canary blind" 's/^\(\[canary\] allocation \)PASS observed=[0-9]*/\1FAIL observed=0/'
    red_check evidence 11 "0 '[verdict] family dsp' lines" "family verdict missing" '/^\[verdict\] family dsp /d'
    red_check evidence 12 "the dsp family did not pass" "family verdict FAIL" '/^\[verdict\] family dsp /s/ verdict=PASS / verdict=FAIL /'
    red_check evidence 13 "names its image kind as 'timing'" "image kind of the other image" 's/^\[meta\] kind=evidence$/[meta] kind=timing/'
    red_check evidence 14 "carries no '[meta] posture' line" "posture line missing" '/^\[meta\] posture /d'
    red_check timing 15 "cycle counter did not prove itself live" "cycle counter not live" 's/^\[canary\] cycle_counter PASS.*/[canary] cycle_counter UNAVAILABLE state=locked/'
    red_check timing 16 "carries no sample count" "timing row without its sample count" 's/^\(\[spread\] dsp\.step,warm,\)[0-9]*,/\1,/'
    red_check timing 17 "carries 2 'dsp.step' cold rows" "timed region with two cold rows" '/^\[spread\] dsp\.step,cold,/p'
    red_check evidence 18 "nonzero harness floor" "watermark floor nonzero" '/^\[verdict\] stack region=startup /s/ floor=0 / floor=8 /'
    red_check evidence 19 "saturated window" "watermark window saturated" '/^\[verdict\] stack region=startup /s/ saturated=no/ saturated=yes/'
    red_check evidence 20 "carries 0 'predictive.steady' stack lines" "stack region missing" '/^\[verdict\] stack region=predictive\.steady /d'
    red_check evidence 21 "exceeds the image's ${reserve} B reserve" "stack high-water past the reserve" "/^\\[verdict\\] stack region=program /s/ from_top=[0-9]* / from_top=$((reserve + 4)) /"
    red_check evidence 22 "prints a reserve of $((reserve - 4)) B" "stack line disagreeing with the image's reserve" "/^\\[verdict\\] stack region=startup /s/ reserve=[0-9]* / reserve=$((reserve - 4)) /"
    red_check evidence 23 "carries no '[spread]' line" "stream tag missing" '/^\[spread\] /d'
    red_check evidence 24 "the field 'x0' is named in both class headers" "field under both class tags" 's/^\(\[spread\] region,.*\)wraps$/\1x0/'
    red_check evidence 25 "carries index 4 where 3 is required" "exact rows reordered" '/^\[exact\] 3,/{h;d};/^\[exact\] 4,/G'
}

self_test_cases_series()
{
    local row="" flipped="" raised=""
    row="$(grep -m 1 '^\[exact\] 100,' "${self_test_root}/evidence.txt")"
    flipped="${row%?}$([ "${row: -1}" = "0" ] && echo 1 || echo 0)"
    row="$(grep -m 1 '^\[spread\] dsp\.step,warm,' "${self_test_root}/timing.txt")"
    raised="$(printf '%s\n' "${row}" | awk -F, -v OFS=, '{ $8 = $8 + 7; print }')"

    red_series evidence pin "" 26 "decoded reset=software where the series requires reset=pin" "run of another reset class" 's/^\(\[meta\] posture reset=\)pin /\1software /'
    red_series evidence pin "" 27 "in the exact-class block at block line" "exact block one hex digit apart" "s|^\\[exact\\] 100,.*|${flipped}|"
    red_series timing pin 6 28 "past the stated bound of 6 cycles" "raw_max 7 cycles apart under a bound of 6" "s|^\\[spread\\] dsp\\.step,warm,.*|${raised}|"
    red_series evidence pin "" 29 "reports the posture" "posture differing outside the reset fields" 's/^\(\[meta\] posture .*\) core_hz=[0-9]*/\1 core_hz=48000000/'

    expect_status 0 "" "raw_max 7 cycles apart under a bound of 7" self_test_compare timing pin 7 "$(series timing "s|^\\[spread\\] dsp\\.step,warm,.*|${raised}|")"
}

# The same refusals on the ESP32 leg's recordings, where its stream differs
# from the NUCLEO-H753ZI's: the scalar's width, the absent family, the posture
# fields, the schedule line, the task stack and the reset classes.
self_test_cases_esp32()
{
    red_check esp32-evidence 6 "is NOT the built one" "ESP32: fingerprint of another image" 's/^\(\[meta\] build_id=\)[0-9a-f]*/\1'"$(printf '0%.0s' $(seq 64))"'/'
    red_check esp32-evidence 11 "must be declared ABSENT on exactly one line; the report carries 0" "ESP32: absent family undeclared" '/^\[verdict\] family predictive /d'
    red_check esp32-evidence 11 "must be declared ABSENT on exactly one line; the report carries 1" "ESP32: absent family claiming a verdict" '/^\[verdict\] family predictive /s/ verdict=ABSENT / verdict=PASS /'
    red_check esp32-evidence 14 "carries no '[meta] schedule' line" "ESP32: schedule line missing" '/^\[meta\] schedule /d'
    red_check esp32-evidence 14 "carries no 'reset_rom=' field" "ESP32: posture without its ROM reset code" 's/ reset_rom=[0-9]*//'
    red_check esp32-timing 17 "carries 0 'trajectory.step' warm rows" "ESP32: timed region without its warm row" '/^\[spread\] trajectory\.step,warm,/d'
    red_check esp32-evidence 19 "used to its last byte" "ESP32: task stack used to its last byte" 's/^\(\[verdict\] stack region=program \).*/\1used=8192 free=0 of=8192 bytes/'
    red_check esp32-evidence 20 "does not add up" "ESP32: task stack figures that do not sum" '/^\[verdict\] stack region=program /s/ used=\([0-9]*\)/ used=\14/'
    red_check esp32-evidence 22 "names a 4096 B stack where the source names 8192 B" "ESP32: task stack disagreeing with the source" 's/^\(\[verdict\] stack region=program \).*/\1used=3152 free=944 of=4096 bytes/'
    red_check esp32-evidence 25 "is not a bit pattern of 8 hexadecimal digits" "ESP32: binary64-width value on the binary32 leg" '/^\[exact\] 7,/s/$/00000000/'
    red_series esp32-evidence power_on "" 26 "decoded reset=panic where the series requires reset=power_on" "ESP32: run of another reset class" 's/^\(\[meta\] posture reset=\)power_on reset_reason=1 /\1panic reset_reason=4 /' power-cycle
    red_series esp32-evidence power_on "" 29 "reports the posture" "ESP32: posture differing outside the reset fields" 's/ core_mhz=160 / core_mhz=240 /' power-cycle
    red_check esp32-evidence 9 "carries no '[canary] heap_view' line" "ESP32: heap-view canary missing" '/^\[canary\] heap_view /d'
    red_check esp32-evidence 9 "carries no self-count" "ESP32: heap-view canary without its self-count" '/^\[canary\] heap_view /s/ self=[0-9]*//'
    red_check esp32-evidence 10 "heap-view canary did not pass" "ESP32: heap-view canary blind" 's/^\(\[canary\] heap_view \)PASS/\1FAIL/'
    red_check esp32-evidence 30 "carries no heap_view verdict" "ESP32: heap figures without a verdict" '/^\[verdict\] family dsp /s/ heap_windows=1 heap_view=clean$//'
    red_check esp32-evidence 30 "heap view is not clean" "ESP32: window kept more than the view's snapshot" '/^\[verdict\] family dsp /s/ heap_view=clean$/ heap_view=FAIL/'

    expect_status 0 "" "ESP32: evidence recording, unchanged" self_test_check esp32-evidence "${self_test_root}/esp32-evidence.txt"
    expect_status 0 "" "ESP32: timing recording, unchanged" self_test_check esp32-timing "${self_test_root}/esp32-timing.txt"
    expect_status 0 "" "ESP32: three unchanged timing runs, reset=power_on, bound 0" self_test_compare esp32-timing power_on 0 "$(series esp32-timing '' power-cycle)"
    expect_status 0 "" "ESP32: report head behind idle-low NUL bytes, cut and checked" self_test_cut esp32-timing
    self_test_cases_en_pin
}

# The EN-pin class decodes as power_on, so what keeps it apart from a power_on
# series is the restart each run's header records.
self_test_cases_en_pin()
{
    local panic='s/^\(\[meta\] posture reset=\)power_on reset_reason=1 /\1panic reset_reason=4 /' dir=""
    expect_status 2 "en_pin is the ESP32 leg's EN-pin restart class" "en_pin on the NUCLEO-H753ZI" self_test_compare timing en_pin 0 "$(series timing '' reset-only)"
    red_series esp32-evidence en_pin "" 26 "decoded reset=panic where the series requires reset=power_on" "ESP32: en_pin run of another reset class" "${panic}" reset-only
    red_series esp32-evidence en_pin "" 26 "records the restart mode 'power-cycle' where this en_pin series requires 'reset-only'" "ESP32: power-cycle run pooled into an en_pin series" 's/^# mode:      reset-only$/# mode:      power-cycle/' reset-only
    red_series esp32-evidence power_on "" 26 "records the restart mode 'reset-only' where this power_on series requires 'power-cycle'" "ESP32: EN-pin run pooled into a power_on series" 's/^# mode:      power-cycle$/# mode:      reset-only/' power-cycle
    red_series esp32-evidence en_pin "" 26 "records the restart mode 'none' where this en_pin series requires 'reset-only'" "ESP32: run without a restart record in an en_pin series" '/^# mode: /d' reset-only

    expect_status 0 "NOT a power" "ESP32: three unchanged timing runs, en_pin, bound 0, summary naming the EN pulse" self_test_compare_summary esp32-timing en_pin 0 "$(series esp32-timing '' reset-only)" "EN-pin pulse through the console bridge's RTS line, NOT a power"
    dir="$(series esp32-evidence '' reset-only)"
    sed -e "${panic}" "${dir}/run-02.txt" > "${dir}/refused-01.txt"
    expect_status 0 "refused-01.txt (refused, not counted) reset=panic" "ESP32: en_pin series listing a run of another class as refused" self_test_compare_summary esp32-evidence en_pin "" "${dir}" "refused-01.txt (refused, not counted) reset=panic"
}

# The boot handshake, checked as a built image's report would be, with nothing
# predated, over the line the NUCLEO-H753ZI evidence recording carries.
self_test_cases_current()
{
    local handshake='s/^\(\[meta\] boot_handshake state=\)received waited_ms=[0-9]*/\1STATE waited_ms=WAITED/'

    red_current evidence 5 "carries no '[meta] boot_handshake' line" "boot handshake line missing" '/^\[meta\] boot_handshake /d'
    red_current evidence 5 "boot handshake was not received" "boot handshake expired" "$(printf '%s' "${handshake}" | sed 's/STATE/timeout/;s/WAITED/30000/')"

    expect_status 0 "" "boot handshake received" self_test_check_current evidence "$(corrupt evidence '')"
}

self_test_cases_green()
{
    expect_status 0 "" "evidence recording, unchanged" self_test_check evidence "${self_test_root}/evidence.txt"
    expect_status 0 "" "timing recording, unchanged" self_test_check timing "${self_test_root}/timing.txt"
    expect_status 0 "" "three unchanged evidence runs, reset=pin" self_test_compare evidence pin "" "$(series evidence '')"
    expect_status 0 "" "three unchanged timing runs, reset=pin, bound 0" self_test_compare timing pin 0 "$(series timing '')"
}

do_self_test()
{
    local file="" status="" implemented="" missing="" counts=""
    for file in evidence.txt timing.txt reference.txt esp32-evidence.txt esp32-timing.txt esp32-reference.txt; do
        [ -f "${fixture_dir}/${file}" ] || fail 3 "the fixture ${fixture_dir}/${file} is missing."
    done
    self_test_root="$(mktemp -d "${TMPDIR:-/tmp}/ctrlpp-board-self-test.XXXXXX")"
    trap 'rm -rf -- "${self_test_root}"' EXIT
    cp -R "${fixture_dir}/." "${self_test_root}/"

    echo "board_capture.sh self-test: committed recordings, replayed from a temporary copy; no board, port, flash or built image is used."
    self_test_cases_green
    self_test_cases_single
    self_test_cases_series
    self_test_cases_esp32
    self_test_cases_current

    implemented="$(implemented_refusals)"
    for status in ${implemented}; do
        printf '%s\n' "${self_test_red[@]+"${self_test_red[@]}"}" | grep -qx "${status}" || missing="${missing} ${status}"
    done
    counts="${#self_test_red[@]} refusal assertions seen red against the $(printf '%s\n' "${implemented}" | wc -l) refusal statuses this script implements ($(printf '%s\n' "${implemented}" | head -1) to $(printf '%s\n' "${implemented}" | tail -1)); ${self_test_green} acceptance assertions seen green"
    if [ "${self_test_failures}" -ne 0 ] || [ -n "${missing}" ]; then
        echo "self-test FAILED: ${counts}; ${self_test_failures} assertion(s) did not hold; statuses with no red case:${missing:- none}." >&2
        exit 1
    fi
    echo "self-test PASS: ${counts}."
}

if [ "${self_test}" = "yes" ]; then
    [ -z "${subcommand}" ] || fail 2 "--self-test replays its own fixtures and takes no subcommand or operand."
    do_self_test
    exit 0
fi

case "${subcommand}" in
    prepare)
        [ -z "${destination}" ] || fail 2 "prepare takes no destination."
        do_prepare
        ;;
    capture)
        do_capture
        ;;
    check)
        do_check
        ;;
    compare)
        do_compare
        ;;
    discover)
        [ -z "${destination}" ] || fail 2 "discover takes no operand."
        do_discover
        ;;
    "")
        usage
        fail 2 "no subcommand."
        ;;
    *)
        fail 2 "unknown subcommand '${subcommand}'. Expected 'prepare', 'capture', 'check', 'compare' or 'discover'."
        ;;
esac
