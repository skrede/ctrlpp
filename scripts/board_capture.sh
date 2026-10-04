#!/usr/bin/env bash
#
# board_capture.sh -- build the NUCLEO-H753ZI image, then capture one boot-once
# report from the board and refuse anything that is not evidence.
#
# Usage:
#   board_capture.sh prepare
#   board_capture.sh capture [--reset-only | --no-reset] <destination-file>
#
# "prepare" configures, builds and produces the raw binary beside the ELF. It
# touches no hardware. "capture" performs no build at all: a build discovered
# inside a capture is how a capture ends up describing an image that was never
# flashed, which is the failure this split exists to prevent, so a missing ELF
# or binary is an error naming the target.
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
# Environment:
#   CAPTURE_BUILD_DIR        build tree (default: <repo>/build-arm-capture)
#   CAPTURE_MOUNT            ST-LINK mass-storage mount
#                            (default: /run/media/$USER/NOD_H753ZI)
#   CAPTURE_PORT             serial port, when discovery must be overridden
#   CAPTURE_READ_TIMEOUT     seconds (default: 45)
#   CAPTURE_BYTE_CAP         bytes (default: 262144)
#   CAPTURE_BUILD_JOBS       build parallelism (default: 6)
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
#   3. A passing capture says the report arrived and matched. It says nothing
#      about how the target was reset, and this script cannot distinguish a power
#      cycle from a warm reset -- that distinction has to come from the firmware's
#      own reset-cause reporting, which this script only transports.
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
#   6. A family's PASS says the value its last step returned lies within its
#      derived bound of the host reference; it does not compare the steps before
#      it. The allocation count on the same line is transported, not judged.

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
target_name="ctrlpp_nucleo_h753zi"
elf_path="${build_dir}/${target_name}.elf"
bin_path="${build_dir}/${target_name}.bin"

# Spelled out because ninja 1.13.2 aborts on this project's dynamic dependency
# file. Parallelism belongs to the build only; the capture is serial.
generator="Unix Makefiles"

# The image's closing line. Seeing it means the report is complete, so the read
# can stop before the timeout rather than always paying it.
report_end_marker="golden diff "

# Every family the image runs. Each must print exactly one verdict line, so a
# family that ran and reported nothing is a refusal rather than a shorter report.
expected_families="control estimation dsp trajectory"

# --- Arguments ----------------------------------------------------------------

subcommand=""
capture_mode="flash"
destination=""

for arg in "$@"; do
    case "${arg}" in
        -h|--help)
            sed -n '2,91p' "${BASH_SOURCE[0]}"
            exit 0
            ;;
        --reset-only)
            capture_mode="reset-only"
            ;;
        --no-reset)
            capture_mode="no-reset"
            ;;
        -*)
            echo "ERROR: unknown option '${arg}'." >&2
            echo "  Usage: scripts/board_capture.sh prepare" >&2
            echo "         scripts/board_capture.sh capture [--reset-only | --no-reset] <destination-file>" >&2
            exit 2
            ;;
        *)
            if [ -z "${subcommand}" ]; then
                subcommand="${arg}"
            elif [ -z "${destination}" ]; then
                destination="${arg}"
            else
                echo "ERROR: unexpected extra argument '${arg}'." >&2
                exit 2
            fi
            ;;
    esac
done

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
        || fail 3 "build of ${target_name} failed."

    [ -f "${elf_path}" ] || fail 3 "the build produced no ${elf_path}."

    "${toolchain_prefix}objcopy" -O binary "${elf_path}" "${bin_path}" \
        || fail 3 "objcopy could not produce ${bin_path}."

    echo "prepared: ${elf_path}"
    echo "prepared: ${bin_path}"
    echo "build id: $(elf_build_id)"
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

elf_build_id()
{
    "${toolchain_prefix}readelf" -n "${elf_path}" \
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

    verify_canary

    grep -q "${report_end_marker}" "${transcript}" \
        || fail 7 "the captured report carries no '${report_end_marker}' verdict line -- it is truncated or the run did not complete."

    if grep -q "${report_end_marker}FAIL" "${transcript}"; then
        fail 8 "the board's golden verdict is FAIL."
    fi

    verify_families

    return 0
}

verify_families()
{
    local family="" line="" count=0
    for family in ${expected_families}; do
        count="$(grep -c "^\[family\] ${family} " "${transcript}" || true)"
        [ "${count}" -eq 1 ] \
            || fail 11 "the captured report carries ${count} '[family] ${family}' lines where exactly one is required."
        line="$(grep -m 1 "^\[family\] ${family} " "${transcript}")"
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
    canary_line="$(grep -m 1 '^\[canary\] ' "${transcript}" || true)"
    [ -n "${canary_line}" ] \
        || fail 9 "the captured report carries no '[canary]' line, so its allocation figure comes from a sensor that never proved it can see an allocation."

    printf '%s\n' "${canary_line}" | grep -q '^\[canary\] PASS observed=[1-9]' \
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
        echo "# build id:  ${built_id}"
        echo "# toolchain: $("${toolchain_prefix}gcc" --version | head -1)"
        echo "#"
        echo "# The build id above is the linker's SHA-1 note, read out of the ELF with"
        echo "# ${toolchain_prefix}readelf -n and matched against the line the firmware printed."
        echo "# It is an identity check against accidental staleness, not a signature."
        echo
        cat "${transcript}"
    } > "${destination}"
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

case "${subcommand}" in
    prepare)
        [ -z "${destination}" ] || fail 2 "prepare takes no destination."
        do_prepare
        ;;
    capture)
        do_capture
        ;;
    "")
        fail 2 "no subcommand. Usage: scripts/board_capture.sh {prepare|capture <destination-file>}"
        ;;
    *)
        fail 2 "unknown subcommand '${subcommand}'. Expected 'prepare' or 'capture'."
        ;;
esac
