#!/usr/bin/env bash
#
# board_capture.sh -- build the NUCLEO-H753ZI image, then capture one boot-once
# report from the board and refuse anything that is not evidence.
#
# Usage:
#   board_capture.sh prepare
#   board_capture.sh capture [--timing] [--reset-only | --no-reset] <destination-file>
#   board_capture.sh check [--timing] <transcript-or-capture-file>
#   board_capture.sh compare [--timing --bound=<cycles>] <count> <reset-class> <destination-directory>
#   board_capture.sh compare --replay [--timing --bound=<cycles>] <reset-class> <series-directory>
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
#                  tool here can cut the target's power. A run whose decoded
#                  class is not the requested one is refused and not counted.
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
#
# The reader is attached BEFORE anything resets the target. The image prints its
# report once at boot and then idles, so a reader attached afterwards reads an
# empty stream; that is a failed capture and this script says so rather than
# reporting a pass it did not earn.
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
#   CAPTURE_MOUNT            ST-LINK mass-storage mount
#                            (default: /run/media/$USER/NOD_H753ZI)
#   CAPTURE_PORT             serial port, when discovery must be overridden
#   CAPTURE_READ_TIMEOUT     seconds (default: 45)
#   CAPTURE_BYTE_CAP         bytes (default: 262144)
#   CAPTURE_BUILD_JOBS       build parallelism (default: 6)
#   COMPARE_ATTEMPTS         captures per run when the report head is lost,
#                            status 4 or 5 (default: 3)
#   COMPARE_PLUG_TIMEOUT     seconds the power-on prompt waits for the node to
#                            vanish, and again to return (default: 120)
#
# Exit status:
#   0  capture written    2  usage       3  environment or build
#   4  empty stream       5  no fingerprint line
#   6  fingerprint mismatch (the resident image is not the built one)
#   7  no verdict line    8  verdict FAIL
#   9  no allocation-sensor canary line
#  10  canary FAIL (the sensor did not observe its deliberate allocations)
#  11  a family's verdict line is missing or repeated
#  12  a family's verdict is not PASS
#  13  the image's kind line is missing or names the other image
#  14  no posture line
#  15  timing UNAVAILABLE: the cycle counter did not prove itself live, or the
#      timing image printed no timing line at all
#  16  a timing line without its sample count, or a per-iteration figure
#      without the block size it was divided by
#  17  a timed region lacks exactly one cold and one warm line
#  18  a stack watermark line reports a nonzero harness floor: the gap is
#      smaller than the harness's own frames, so its figure is not a
#      measurement of the chain
#  19  a stack watermark line reports a saturated window: the chain reached the
#      window's bottom word, so its depth is unknown rather than small
#  20  a stack watermark line lacks its floor, its saturation verdict, its
#      figure or its reserve, or a measured region's line is missing or repeated
#  21  a stack line's high-water exceeds the image's stack reserve: the
#      measured stack needs more than the linker reserved for it
#  22  the reserve a stack line prints disagrees with the _Min_Stack_Size
#      symbol in the ELF, so the line and the image disagree on the reserve
#  23  one of the five stream tags is missing, so a block of the report is
#      absent rather than empty
#  24  a field name appears in both the exact-class and the spread-class
#      header, or twice in one, so it would be claimed bitwise identical and
#      bounded by a spread at once
#  25  a class block is malformed: its header does not key its rows (the loop
#      index for the exact class, the region for the spread class), a row's
#      width differs from its header's, an exact row's index breaks the
#      sequence 0, 1, 2, ... (a reordering is a difference, not a pass), or an
#      exact value is not a 16-digit bit pattern
#  26  a run's decoded reset class is not the class the series requires
#  27  the exact-class blocks of two runs differ
#  28  a cycle field's spread across the runs exceeds the stated bound
#  29  the runs do not share what a comparison needs them to share: the build
#      id, the posture outside the reset fields, the spread rows and their
#      shape fields, or a spread column this script can classify
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
#      The reader attaches after the ST-LINK enumerates and the target does not
#      wait for that, so a power-on report can lose its head and refuse with
#      status 5.
#
#  11. The spread bound is one absolute cycle count applied to every timed
#      region's raw figures, from the 129-cycle overhead to the predictive
#      solve; it is the invocation's stated claim, not a figure derived here.

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

# --- Arguments ----------------------------------------------------------------

subcommand=""
capture_mode="flash"
image_kind="evidence"
destination=""
replay="no"
spread_bound=""
operands=()

usage()
{
    echo "  Usage: scripts/board_capture.sh prepare" >&2
    echo "         scripts/board_capture.sh capture [--timing] [--reset-only | --no-reset] <destination-file>" >&2
    echo "         scripts/board_capture.sh check [--timing] <transcript-or-capture-file>" >&2
    echo "         scripts/board_capture.sh compare [--timing --bound=<cycles>] <count> <reset-class> <destination-directory>" >&2
    echo "         scripts/board_capture.sh compare --replay [--timing --bound=<cycles>] <reset-class> <series-directory>" >&2
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
        --replay)
            replay="yes"
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
    cmake -S "${leg_dir}" -B "${build_dir}" \
        -G "${generator}" \
        -DCMAKE_TOOLCHAIN_FILE="${leg_dir}/cmake/toolchain-arm-none-eabi.cmake" \
        -DCMAKE_BUILD_TYPE=Release \
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

# --- discovery ----------------------------------------------------------------

resolve_port()
{
    if [ -n "${CAPTURE_PORT:-}" ]; then
        printf '%s\n' "${CAPTURE_PORT}"
        return 0
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
    value="$("${toolchain_prefix}nm" "${elf_path}" | awk '$3 == "_Min_Stack_Size" { print $1 }')"
    [ -n "${value}" ] || return 0
    printf '%d\n' "0x${value}"
}

elf_build_id()
{
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
    if [ -n "${reader_pid}" ]; then
        kill "${reader_pid}" 2>/dev/null || true
        wait "${reader_pid}" 2>/dev/null || true
        reader_pid=""
    fi
    if [ -n "${raw_stream}" ]; then
        rm -f -- "${raw_stream}"
    fi
    if [ -n "${transcript}" ]; then
        rm -f -- "${transcript}"
    fi
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
    local waited=0 size=0
    while [ "${waited}" -lt "${read_timeout}" ]; do
        size="$(wc -c < "${raw_stream}" 2>/dev/null || echo 0)"
        if [ "${size}" -ge "${byte_cap}" ]; then
            return 0
        fi
        if grep -q "${report_end_marker}" "${raw_stream}" 2>/dev/null; then
            return 0
        fi
        sleep 1
        waited=$((waited + 1))
    done
    return 0
}

# --- refusals -----------------------------------------------------------------

verify_transcript()
{
    local built_id="$1" printed_id=""

    grep -q '[^[:space:]]' "${transcript}" \
        || fail 4 "the capture read an EMPTY stream from the board -- no report was observed, and an absent report is not a pass. The image reports once at boot: reset the target while this script holds the port (drop the --no-reset flag), or press B2 and re-run."

    printed_id="$(sed -n 's/^\[meta\] build_id=\([0-9a-f]*\).*$/\1/p' "${transcript}" | head -1)"
    [ -n "${printed_id}" ] \
        || fail 5 "the captured report carries no '[meta] build_id=' line, so the image that produced it cannot be identified."

    [ "${printed_id}" = "${built_id}" ] \
        || fail 6 "the running image is NOT the built one: the board reported ${printed_id} but ${elf_path} is ${built_id}. Flash the target before capturing."

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
    verify_tags
    verify_spread_header

    if [ "${image_kind}" = "timing" ]; then
        verify_timing
    else
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
# bitwise claim.
verify_exact()
{
    local problem=""
    problem="$(block_lines exact | awk -F, '
        NR == 1 { width = NF; if ($1 != "k") { print "the exact-class header does not key its rows by the loop index: \047" $0 "\047"; exit } next }
        NF != width { print "an exact-class row carries " NF " fields against a header of " width ": \047" $0 "\047"; exit }
        $1 != (NR - 2) "" { print "an exact-class row carries index " $1 " where " NR - 2 " is required: \047" $0 "\047"; exit }
        { for (i = 2; i <= NF; ++i) if ($i !~ /^[0-9A-F]{16}$/) { print "an exact-class value is not a 16-digit bit pattern: \047" $0 "\047"; exit } }
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
    for field in reset reset_flags reset_rsr core_hz icache dcache flash_latency_ws cycle_wrap_ms compiler; do
        printf '%s\n' "${posture}" | grep -Eq " ${field}=[^ ]+" \
            || fail 14 "the posture line carries no '${field}=' field: '${posture}'."
    done
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
    reserve="$(elf_stack_reserve)"
    [ -n "${reserve}" ] || fail 3 "${elf_path} carries no _Min_Stack_Size symbol, so the stack reserve is unknown."
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
            || fail 22 "a stack line prints a reserve of ${printed} B where ${elf_path} has _Min_Stack_Size = ${reserve} B: '${line}'."
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

# --- capture ------------------------------------------------------------------

do_capture()
{
    [ -n "${destination}" ] || fail 2 "capture needs a destination file."
    [ -d "$(dirname "${destination}")" ] \
        || fail 2 "the destination directory '$(dirname "${destination}")' does not exist."
    [ -f "${elf_path}" ] || fail 3 "no ${elf_path}; run 'board_capture.sh prepare' first."
    [ -f "${bin_path}" ] || fail 3 "no ${bin_path}; run 'board_capture.sh prepare' first."

    local built_id port
    built_id="$(elf_build_id)"
    [ -n "${built_id}" ] || fail 3 "${elf_path} carries no NT_GNU_BUILD_ID note; the image cannot be fingerprinted."

    port="$(resolve_port)" || fail 3 "no ST-LINK serial port found: neither /dev/serial/by-id/*STLINK* nor a USB device 0483:374e with a tty child resolved."
    [ -c "${port}" ] || fail 3 "'${port}' is not a character device."

    if [ "${capture_mode}" = "flash" ]; then
        [ -d "${mount_dir}" ] || fail 3 "'${mount_dir}' is not a directory; the ST-LINK mass storage is not mounted."
        mountpoint -q "${mount_dir}" || fail 3 "'${mount_dir}' exists but is not a mount point; refusing to write the image into the local filesystem."
    fi

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

    wait_for_report
    reader_stop

    tr -d '\r' < "${raw_stream}" > "${transcript}"
    verify_transcript "${built_id}"
    write_artifact "${port}" "${built_id}"

    echo "capture written: ${destination}"
    echo "build id:        ${built_id}"
}

# The same refusals over a file already on disk. A capture file's header lines
# start with '#', which no refusal's anchored pattern matches.
do_check()
{
    [ -n "${destination}" ] || fail 2 "check needs a transcript file."
    [ -f "${destination}" ] || fail 2 "'${destination}' is not a file."
    [ -f "${elf_path}" ] || fail 3 "no ${elf_path}; run 'board_capture.sh prepare' first."

    local built_id
    built_id="$(elf_build_id)"
    transcript="$(mktemp "${TMPDIR:-/tmp}/ctrlpp-board-check.XXXXXX")"
    trap cleanup EXIT
    tr -d '\r' < "${destination}" > "${transcript}"
    verify_transcript "${built_id}"
    echo "check passed: ${destination}"
}

# --- compare ------------------------------------------------------------------
#
# A series is N restarts of one resident image, each captured and refused on its
# own, then compared: the exact-class block byte for byte, the spread-class block
# field by field against a stated bound. The capture subcommand does the
# capturing, so every run of a series passes every refusal a single capture does.

script_path="${script_dir}/$(basename "${BASH_SOURCE[0]}")"
compare_attempts="${COMPARE_ATTEMPTS:-3}"
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
}

# A lost report head (an empty stream, or no fingerprint line) is the transport
# failing, not the run, so it is retried; any other refusal ends the series.
collect_one()
{
    local file="$1" mode="$2" attempt=1 status=0
    while [ "${attempt}" -le "${compare_attempts}" ]; do
        status=0
        if [ "${mode}" = "power_on" ]; then
            await_power_cycle
            "${script_path}" capture ${kind_flag} --no-reset "${file}" || status=$?
        else
            "${script_path}" capture ${kind_flag} --reset-only "${file}" || status=$?
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
collect_series()
{
    local index=1 file="" status=0 class=""
    collect_one "${series_dir}/priming.txt" pin || fail "$?" "the priming run was refused; the series was not started."
    while [ "${index}" -le "${series_count}" ]; do
        file="$(printf '%s/run-%02d.txt' "${series_dir}" "${index}")"
        status=0
        collect_one "${file}" "${requested_class}" || status=$?
        [ "${status}" -eq 0 ] || fail "${status}" "run ${index} of the series was refused by the capture; the series stops at $((index - 1)) runs."
        class="$(posture_field "${file}" reset)"
        if [ "${class}" != "${requested_class}" ]; then
            mv -- "${file}" "${file%.txt}-refused.txt"
            compare_fail 26 "run ${index} decoded reset=${class:-absent} where the series requires reset=${requested_class}; it is kept as $(basename "${file%.txt}")-refused.txt and not counted, and the series stops."
            return 0
        fi
        index=$((index + 1))
    done
    return 0
}

check_classes()
{
    local run="" class=""
    for run in "${runs[@]}"; do
        class="$(posture_field "${run}" reset)"
        [ "${class}" = "${requested_class}" ] \
            || compare_fail 26 "$(basename "${run}") decoded reset=${class:-absent} where the series requires reset=${requested_class}."
    done
    return 0
}

# The fingerprint and every posture field but the reset's own must agree, since
# a spread is only reported against the one posture that produced it.
check_shared()
{
    local run="" first_id="" id="" first_posture="" posture=""
    first_id="$(block_of meta "${runs[0]}" | sed -n 's/^\[meta\] build_id=\([0-9a-f]*\).*/\1/p' | head -1)"
    first_posture="$(posture_of "${runs[0]}" | sed -E 's/ reset(_flags|_rsr)?=[^ ]*//g')"
    for run in "${runs[@]}"; do
        id="$(block_of meta "${run}" | sed -n 's/^\[meta\] build_id=\([0-9a-f]*\).*/\1/p' | head -1)"
        [ "${id}" = "${first_id}" ] \
            || compare_fail 29 "$(basename "${run}") carries build id ${id:-absent} where $(basename "${runs[0]}") carries ${first_id:-absent}; a series is one image."
        posture="$(posture_of "${run}" | sed -E 's/ reset(_flags|_rsr)?=[^ ]*//g')"
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
    local run="" verdict="PASS"
    [ "${compare_status}" -eq 0 ] || verdict="FAIL (status ${compare_status}): ${compare_note}"
    {
        echo "# ctrlpp NUCLEO-H753ZI determinism comparison"
        echo "# compared:      $(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo "# mode:          $([ "${replay}" = "yes" ] && echo "replay of transcripts on disk" || echo "collected by this invocation")"
        echo "# kind:          ${image_kind}"
        echo "# reset class:   ${requested_class} (requested; every counted run's decoded class is listed below)"
        echo "# runs:          ${#runs[@]} counted$([ -n "${series_count}" ] && echo " of ${series_count} requested")"
        echo "# build id:      ${shared_id:-absent} (shared by every counted run)"
        if [ "${image_kind}" = "timing" ]; then
            echo "# evidence:      $(elf_build_id "${evidence_elf}") (the allocation-evidence image built from the same sources)"
        fi
        echo "# toolchain:     $("${toolchain_prefix}gcc" --version | head -1)"
        echo "# posture:       ${shared_posture:-absent}"
        echo "# host load:     $(cut -d' ' -f1-3 /proc/loadavg 2>/dev/null || echo unknown) (the figures count core cycles on the target, not host time)"
        echo "#"
        echo "# Per run, as each image decoded its own reset status register:"
        if [ -f "${series_dir}/priming.txt" ] && [ "${replay}" = "no" ]; then
            echo "#   priming.txt (not counted) $(posture_of "${series_dir}/priming.txt" | grep -o 'reset=[^ ]* reset_flags=[^ ]* reset_rsr=[^ ]*')"
        fi
        for run in "${runs[@]}"; do
            echo "#   $(basename "${run}") $(posture_of "${run}" | grep -o 'reset=[^ ]* reset_flags=[^ ]* reset_rsr=[^ ]*')"
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

compare_arguments()
{
    if [ "${replay}" = "yes" ]; then
        [ "${#operands[@]}" -eq 2 ] || { usage; fail 2 "compare --replay takes a reset class and a series directory."; }
        requested_class="${operands[0]}"
        series_dir="${operands[1]}"
    else
        [ "${#operands[@]}" -eq 3 ] || { usage; fail 2 "compare takes a run count, a reset class and a destination directory."; }
        series_count="${operands[0]}"
        requested_class="${operands[1]}"
        series_dir="${operands[2]}"
        printf '%s\n' "${series_count}" | grep -Eqx '[1-9][0-9]*' && [ "${series_count}" -ge 2 ] \
            || fail 2 "the run count must be an integer of at least 2, since one run compares with nothing: '${series_count}'."
        case "${requested_class}" in
            pin|power_on) ;;
            *) fail 2 "this station can collect a 'pin' or a 'power_on' series, not '${requested_class}'; a series of another class can be compared with --replay." ;;
        esac
        printf '%s\n' "${plug_timeout}" | grep -Eqx '[1-9][0-9]*' \
            || fail 2 "COMPARE_PLUG_TIMEOUT must be a whole number of seconds: '${plug_timeout}'."
        printf '%s\n' "${compare_attempts}" | grep -Eqx '[1-9][0-9]*' \
            || fail 2 "COMPARE_ATTEMPTS must be a positive integer: '${compare_attempts}'."
    fi
    printf '%s\n' "${requested_class}" | grep -Eqx '[a-z_]+' || fail 2 "'${requested_class}' is not a reset class name."
    if [ "${image_kind}" = "timing" ]; then
        printf '%s\n' "${spread_bound}" | grep -Eqx '[0-9]+' \
            || fail 2 "a timing series needs its spread bound stated as --bound=<cycles>, a non-negative integer: '${spread_bound}'."
    else
        [ -z "${spread_bound}" ] || fail 2 "--bound applies to the timing image only; the evidence image times nothing."
    fi
    kind_flag=""
    [ "${image_kind}" = "timing" ] && kind_flag="--timing"
    return 0
}

do_compare()
{
    compare_arguments
    [ -f "${elf_path}" ] || fail 3 "no ${elf_path}; run 'board_capture.sh prepare' first."
    compare_tmp="$(mktemp -d "${TMPDIR:-/tmp}/ctrlpp-board-compare.XXXXXX")"
    trap 'cleanup; rm -rf -- "${compare_tmp}"' EXIT

    local run=""
    if [ "${replay}" = "yes" ]; then
        [ -d "${series_dir}" ] || fail 2 "'${series_dir}' is not a directory."
        for run in "${series_dir}"/run-[0-9][0-9].txt; do
            [ -f "${run}" ] || continue
            "${script_path}" check ${kind_flag} "${run}" > /dev/null || fail "$?" "$(basename "${run}") fails a single-capture refusal; it is not compared."
            runs+=("${run}")
        done
        [ "${#runs[@]}" -ge 2 ] || fail 2 "'${series_dir}' holds ${#runs[@]} run-NN.txt files; a comparison needs at least 2."
    else
        mkdir -p -- "${series_dir}" || fail 2 "cannot create '${series_dir}'."
        ! compgen -G "${series_dir}/run-*.txt" > /dev/null \
            || fail 2 "'${series_dir}' already holds run files; a series is collected into a directory of its own."
        collect_series
        for run in "${series_dir}"/run-[0-9][0-9].txt; do
            [ -f "${run}" ] && runs+=("${run}")
        done
    fi

    check_classes
    if [ "${#runs[@]}" -gt 0 ]; then
        check_shared
        compare_exact
        compare_spread
    fi
    write_summary
    echo "summary written: ${series_dir}/compare-summary.txt"
    [ "${compare_status}" -eq 0 ] || exit "${compare_status}"
    echo "series PASS: ${#runs[@]} runs, reset=${requested_class}, build id ${shared_id}"
}

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
    "")
        usage
        fail 2 "no subcommand."
        ;;
    *)
        fail 2 "unknown subcommand '${subcommand}'. Expected 'prepare', 'capture', 'check' or 'compare'."
        ;;
esac
