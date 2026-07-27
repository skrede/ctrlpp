#!/usr/bin/env bash
#
# Local fuzz smoke run over the built fuzz targets, mirroring the loop the fuzz
# workflow runs so the curated corpus is verifiable without pushing anything.
#
# Each target gets up to two phases, reported separately so a regression is
# never confused with a new discovery:
#
#   Replay   A zero-run pass over the curated seed corpus in
#            tests/fuzz/corpus/<target>, for the targets that have one. This is
#            the deterministic gate: every curated seed must replay without a
#            finding. It is a zero-run pass on purpose, because handing a corpus
#            directory to a fuzzing run makes libFuzzer write newly discovered
#            inputs into it, which would turn a versioned directory into a
#            corpus dump and silently change the regression set.
#
#   Explore  The exploratory run. It is intentionally seedless and short: a
#            smoke test that the target still runs, not a campaign, and it is
#            never handed the curated directory.
#
# Findings are written under the build tree rather than the working tree, so a
# crash artifact never appears in a git status.
#
# Usage:
#   scripts/fuzz_smoke.sh [build-dir]
#   FUZZ_BUILD_DIR=build-fuzz scripts/fuzz_smoke.sh
#
# Environment:
#   FUZZ_BUILD_DIR         build directory holding tests/fuzz (default: build)
#   FUZZ_MAX_TOTAL_TIME    seconds per exploratory run (default: 5)
#   FUZZ_MAX_LEN           input length cap for the exploratory run (default: 1024)

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
cd "${repo_root}"

build_dir="${1:-${FUZZ_BUILD_DIR:-build}}"
max_total_time="${FUZZ_MAX_TOTAL_TIME:-5}"
max_len="${FUZZ_MAX_LEN:-1024}"

if [ ! -d "${build_dir}" ]; then
    echo "ERROR: build directory '${build_dir}' does not exist." >&2
    echo "  Configure and build the fuzz targets first, for example:" >&2
    echo "    cmake -S . -B ${build_dir} -G \"Unix Makefiles\" -DCMAKE_BUILD_TYPE=Debug \\" >&2
    echo "      -DCTRLPP_BUILD_FUZZ_TESTS=ON -DCTRLPP_CMAKE_FETCH_DEPS=ON -DCMAKE_CXX_COMPILER=clang++" >&2
    echo "    cmake --build ${build_dir}" >&2
    exit 1
fi

build_abs="$(cd "${build_dir}" && pwd)"
fuzz_bin_dir="${build_abs}/tests/fuzz"
corpus_root="${repo_root}/tests/fuzz/corpus"
artifact_dir="${build_abs}/fuzz-artifacts"
log_file="${build_abs}/fuzz-smoke.log"

if [ ! -d "${fuzz_bin_dir}" ]; then
    echo "ERROR: no fuzz binaries under '${fuzz_bin_dir}'." >&2
    echo "  Reconfigure with -DCTRLPP_BUILD_FUZZ_TESTS=ON and build." >&2
    exit 1
fi

rm -rf "${artifact_dir}"
mkdir -p "${artifact_dir}"

# An outer wall-clock guard, when the platform provides one. libFuzzer already
# bounds itself, so a missing timeout(1) is not fatal.
timeout_cmd=""
if command -v timeout >/dev/null 2>&1; then
    timeout_cmd="timeout"
fi

# Runs a fuzz binary under the wall-clock guard when available.
run_fuzzer()
{
    local wall_limit="$1"
    shift
    if [ -n "${timeout_cmd}" ]; then
        "${timeout_cmd}" "${wall_limit}" "$@"
    else
        "$@"
    fi
}

# Reports whether a phase left a crash, out-of-memory, timeout, or leak artifact.
artifacts_found()
{
    [ -n "$(ls -A "$1" 2>/dev/null)" ]
}

echo "Build directory:  ${build_abs}"
echo "Corpus root:      ${corpus_root}"
echo "Artifacts:        ${artifact_dir}"
echo

targets=0
failed=0

for fuzz in "${fuzz_bin_dir}"/fuzz_*
do
    [ -x "${fuzz}" ] || continue
    [ -f "${fuzz}" ] || continue

    name="$(basename "${fuzz}")"
    targets=$((targets + 1))
    echo "=== ${name} ==="

    corpus="${corpus_root}/${name}"
    if [ -d "${corpus}" ]; then
        seeds="$(find "${corpus}" -type f | wc -l | tr -d '[:space:]')"
        replay_artifacts="${artifact_dir}/${name}-replay"
        mkdir -p "${replay_artifacts}"
        if run_fuzzer 120 "${fuzz}" -runs=0 \
            "-artifact_prefix=${replay_artifacts}/" "${corpus}" > "${log_file}" 2>&1 \
            && ! artifacts_found "${replay_artifacts}"
        then
            echo "  replay:  PASS (${seeds} curated seed(s))"
            rmdir "${replay_artifacts}" 2>/dev/null || true
        else
            echo "  replay:  FAIL (${seeds} curated seed(s))"
            tail -20 "${log_file}"
            failed=$((failed + 1))
        fi
    fi

    # Seedless by design: the curated directory is deliberately not passed here.
    explore_artifacts="${artifact_dir}/${name}-explore"
    mkdir -p "${explore_artifacts}"
    if run_fuzzer $((max_total_time * 2 + 5)) "${fuzz}" \
        "-max_total_time=${max_total_time}" "-max_len=${max_len}" \
        "-artifact_prefix=${explore_artifacts}/" > "${log_file}" 2>&1 \
        && ! artifacts_found "${explore_artifacts}"
    then
        echo "  explore: PASS (${max_total_time}s, seedless)"
        rmdir "${explore_artifacts}" 2>/dev/null || true
    else
        echo "  explore: FAIL (${max_total_time}s, seedless)"
        tail -20 "${log_file}"
        failed=$((failed + 1))
    fi
    echo
done

rm -f "${log_file}"

if [ "${targets}" -eq 0 ]; then
    echo "ERROR: no fuzz target binaries found under ${fuzz_bin_dir}." >&2
    exit 1
fi

if [ "${failed}" -gt 0 ]; then
    echo "${failed} phase(s) failed across ${targets} target(s). Artifacts under ${artifact_dir}."
    exit 1
fi

rm -rf "${artifact_dir}"
echo "All phases PASS across ${targets} target(s)."
