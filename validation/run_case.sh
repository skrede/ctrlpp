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
if [ -f "${CASE_DIR}/tolerance.cfg" ]; then
    # shellcheck source=/dev/null
    source "${CASE_DIR}/tolerance.cfg"
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
    "$REFERENCE_CSV" "$CANDIDATE_CSV" "$atol" "$rtol" "$ANALYSIS_DIR"
