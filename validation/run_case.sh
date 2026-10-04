#!/usr/bin/env bash
# run_case.sh -- Run one validation case and compare it against its reference.
#
# Usage: ./run_case.sh <case_name> <case_executable>
#   case_name        : directory name under cases/
#   case_executable  : path to the built candidate binary for that case
#
# Writes the reference CSV, the candidate CSV, the plots and the report into
# cases/<case_name>/analysis/, and exits with the comparison's own status.
#
# A case's optional tolerance.cfg is sourced as shell. It may set the case-wide
# pair atol and rtol, and per-column overrides atol_<column> and rtol_<column>,
# where <column> is a CSV header; a column without an override keeps the
# case-wide pair.
#
# Environment:
#   CTRLPP_VALIDATE_OCTAVE  interpreter to invoke (default: octave)
#
# Prerequisites:
#   - the interpreter named above, with the control, signal, splines and
#     quaternion packages
#   - the case executable, already built

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

if [ $# -ne 2 ]; then
    echo "usage: run_case.sh <case_name> <case_executable>" >&2
    exit 2
fi

CASE_NAME="$1"
CASE_BINARY="$2"
OCTAVE="${CTRLPP_VALIDATE_OCTAVE:-octave}"

CASE_DIR="${SCRIPT_DIR}/cases/${CASE_NAME}"
REFERENCE_SCRIPT="${CASE_DIR}/${CASE_NAME}.m"

if [ ! -f "$REFERENCE_SCRIPT" ]; then
    echo "${CASE_NAME}: no reference script at ${REFERENCE_SCRIPT}" >&2
    exit 1
fi

if [ ! -x "$CASE_BINARY" ]; then
    echo "${CASE_NAME}: case executable missing or not executable: ${CASE_BINARY}" >&2
    exit 1
fi

ANALYSIS_DIR="${CASE_DIR}/analysis"
mkdir -p "$ANALYSIS_DIR"

REFERENCE_CSV="${ANALYSIS_DIR}/${CASE_NAME}_octave.csv"
CANDIDATE_CSV="${ANALYSIS_DIR}/${CASE_NAME}_cpp.csv"

atol="1e-10"
rtol="1e-8"
overrides=()
TOLERANCE_CFG="${CASE_DIR}/tolerance.cfg"
if [ -f "$TOLERANCE_CFG" ]; then
    # shellcheck source=/dev/null
    source "$TOLERANCE_CFG"
    # Read in an empty environment, so an atol_* or rtol_* the caller happens to
    # export cannot pass for one the case set.
    while IFS= read -r override; do
        overrides+=("$override")
    done < <(env -i "$BASH" --norc --noprofile -c '
        source "$1"
        for name in $(compgen -v); do
            case "$name" in
                atol_?*|rtol_?*) printf "%s=%s\n" "$name" "${!name}" ;;
            esac
        done' run_case "$TOLERANCE_CFG")
fi

if ! "$OCTAVE" --no-gui "$REFERENCE_SCRIPT" > "$REFERENCE_CSV"; then
    echo "${CASE_NAME}: the reference script failed" >&2
    exit 1
fi

if ! "$CASE_BINARY" > "$CANDIDATE_CSV"; then
    echo "${CASE_NAME}: the case executable failed" >&2
    exit 1
fi

"$OCTAVE" --no-gui "${SCRIPT_DIR}/validate_compare.m" \
    "$REFERENCE_CSV" "$CANDIDATE_CSV" "$atol" "$rtol" "$ANALYSIS_DIR" ${overrides[@]+"${overrides[@]}"}
