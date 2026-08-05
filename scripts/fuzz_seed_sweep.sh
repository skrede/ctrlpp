#!/usr/bin/env bash
#
# Multi-seed and seedless fuzz campaign over the built fuzz targets.
#
# This driver sits beside scripts/fuzz_smoke.sh rather than replacing it. The
# two answer different questions:
#
#   fuzz_smoke.sh   A pre-push mirror of the verdict the fuzz workflow reaches.
#                   It replays the curated corpus deterministically and then
#                   explores briefly, so a regression is caught before a push.
#
#   this script     The evidence a posture claim rests on. A default-seed pass
#                   establishes only that one path through the search is clean,
#                   so every target is swept across a written-down set of
#                   explicit seeds AND run several times with no seed at all.
#                   The seeded runs are reproducible and the seedless ones are
#                   not by construction, which is the point: they explore a
#                   different part of the search than any fixed seed does.
#
# COST IS THE EXECUTION COUNT, NOT THE WALL TIME. Every run is bounded by a
# fixed number of executions rather than by a time budget. A fixed count is
# deterministic and comparable across runs, machines and build configurations;
# a time budget measures how loaded the station was, and an executions-per-
# second figure is a reading of the same contention. Wall time is reported
# separately and labeled, so it stays available for capacity planning without
# being mistaken for a cost.
#
# There is deliberately no outer wall-clock guard. One would silently truncate
# the execution count on a loaded station and take the cost measure with it.
# The bound on a stuck input is libFuzzer's own per-input timeout, whose expiry
# leaves an artifact and is reported as a finding like any other.
#
# Targets listed in tests/fuzz/known_red.txt carry an unsound oracle and a
# named repair. That list exists once and is read by both of its consumers --
# the fuzz workflow and scripts/fuzz_smoke.sh -- and this script reads the same
# one rather than carrying a second copy. A listed target is explored and
# reported like any other: a finding there is a KNOWN RED rather than a
# failure, and a listed target that produces none across the whole campaign is
# reported as a STALE WAIVER.
#
# The curated corpus directory is never handed to a run here. Handing a
# versioned directory to a fuzzing run makes libFuzzer write newly discovered
# inputs into it, which turns the regression set into a corpus dump and changes
# it silently. Corpus replay is the smoke script's business; this one explores.
#
# Findings are kept rather than cleaned up, and each is reported with the path
# to the input that reproduces it. They are written under the build tree rather
# than the working tree, so a crash artifact never appears in a git status.
#
# Usage:
#   scripts/fuzz_seed_sweep.sh [build-dir] [target ...]
#
# Naming targets restricts the sweep to them and leaves every other target's
# output from an earlier invocation alone, which is what makes a single target
# re-triageable at a larger execution count without re-measuring the rest. A
# restricted invocation reports only the targets it swept, and its summary says
# how many those were, so a partial run cannot read as a whole-tree posture.
#
# Environment:
#   FUZZ_BUILD_DIR             build directory holding tests/fuzz (default: build)
#   FUZZ_SWEEP_RUNS            executions per run (default: 100000)
#   FUZZ_SWEEP_SEEDS           explicit seeds, whitespace separated
#                              (default: 1 2 3 4 5 6 7 8)
#   FUZZ_SWEEP_SEEDLESS_RUNS   runs with no seed supplied at all (default: 4)
#   FUZZ_SWEEP_MAX_LEN         input length cap (default: 1024)
#   FUZZ_SWEEP_JOBS            targets swept concurrently (default: 1).
#                              Concurrency cannot move an execution count, so it
#                              cannot move the cost this script measures. It
#                              moves only the wall time, which is already
#                              reported as not comparable across runs.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
cd "${repo_root}"

build_dir="${1:-${FUZZ_BUILD_DIR:-build}}"
[ "$#" -gt 0 ] && shift
requested_targets=("$@")
runs="${FUZZ_SWEEP_RUNS:-100000}"
seeds="${FUZZ_SWEEP_SEEDS:-1 2 3 4 5 6 7 8}"
seedless_runs="${FUZZ_SWEEP_SEEDLESS_RUNS:-4}"
max_len="${FUZZ_SWEEP_MAX_LEN:-1024}"
jobs="${FUZZ_SWEEP_JOBS:-1}"

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
known_red_file="${repo_root}/tests/fuzz/known_red.txt"
campaign_dir="${build_abs}/fuzz-campaign"

if [ ! -d "${fuzz_bin_dir}" ]; then
    echo "ERROR: no fuzz binaries under '${fuzz_bin_dir}'." >&2
    echo "  Reconfigure with -DCTRLPP_BUILD_FUZZ_TESTS=ON and build." >&2
    exit 1
fi

mkdir -p "${campaign_dir}"

# The reason a target is waived, empty for a target that is not. Read from the
# same list the fuzz workflow and the smoke script read, so no two consumers can
# disagree about which targets are waived.
known_red_target()
{
    [ -f "${known_red_file}" ] || return 0
    awk -v target="$1" '
      /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
      $1 == target { $1 = ""; sub(/^[[:space:]]+/, ""); print; exit }
    ' "${known_red_file}"
}

# Whole seconds of wall time between two nanosecond stamps, for the reported
# figure only. It is never used to bound a run.
elapsed_seconds()
{
    echo $((($2 - $1) / 1000000000))
}

# One leg of one target: a single bounded run, its artifacts kept under a
# directory of its own so a finding is attributable to the run that produced it.
# Prints the run's report lines and returns the number of findings it left.
run_leg()
{
    local fuzz="$1"
    local label="$2"
    local leg_dir="$3"
    shift 3

    mkdir -p "${leg_dir}"
    local started
    started="$(date +%s%N)"
    local status=0
    "${fuzz}" "-runs=${runs}" "-max_len=${max_len}" "$@" \
        "-artifact_prefix=${leg_dir}/" > "${leg_dir}/output.log" 2>&1 || status=$?
    local finished
    finished="$(date +%s%N)"

    local found=0
    local artifact
    for artifact in "${leg_dir}"/*
    do
        [ -f "${artifact}" ] || continue
        case "${artifact}" in
            */output.log) continue ;;
        esac
        echo "    ${label}: FINDING -- reproduce with"
        echo "      ${fuzz} ${artifact}"
        found=$((found + 1))
    done

    if [ "${found}" -eq 0 ] && [ "${status}" -ne 0 ]; then
        # A nonzero exit with no artifact is a defect in the run itself rather
        # than a discovery, and it is reported rather than counted as clean.
        echo "    ${label}: EXITED ${status} WITHOUT AN ARTIFACT -- see ${leg_dir}/output.log"
        found=$((found + 1))
    fi

    echo "$(elapsed_seconds "${started}" "${finished}") ${label}" >> "${leg_dir}/../wall.txt"
    return "${found}"
}

# The whole campaign for one target, written to a report file so concurrent
# workers cannot interleave their output and the final report keeps a fixed
# target order.
sweep_target()
{
    local name="$1"
    local fuzz="${fuzz_bin_dir}/${name}"
    local target_dir="${campaign_dir}/${name}"
    local report="${campaign_dir}/${name}.report"
    local waiver
    waiver="$(known_red_target "${name}")"

    # Only this target's output is discarded, so naming a subset re-measures the
    # subset and leaves the rest of an earlier invocation intact.
    rm -rf "${target_dir}" "${report}" "${campaign_dir}/${name}.status"
    mkdir -p "${target_dir}"
    : > "${target_dir}/wall.txt"

    local findings=0
    local seeded_count=0
    local seedless_count=0
    local leg_findings

    {
        if [ -n "${waiver}" ]; then
            echo "=== ${name} === WAIVED -- ${waiver}"
        else
            echo "=== ${name} ==="
        fi

        local seed
        for seed in ${seeds}
        do
            seeded_count=$((seeded_count + 1))
            leg_findings=0
            run_leg "${fuzz}" "seed ${seed}" "${target_dir}/seed-${seed}" "-seed=${seed}" || leg_findings=$?
            findings=$((findings + leg_findings))
        done

        local index
        for index in $(seq 1 "${seedless_runs}")
        do
            seedless_count=$((seedless_count + 1))
            leg_findings=0
            run_leg "${fuzz}" "seedless run ${index}" "${target_dir}/seedless-${index}" || leg_findings=$?
            findings=$((findings + leg_findings))
        done

        echo "  seeded:    ${seeded_count} run(s) at -runs=${runs} each, seeds ${seeds}"
        echo "  seedless:  ${seedless_count} run(s) at -runs=${runs} each, no seed supplied"
        echo "  executions: $(( (seeded_count + seedless_count) * runs )) total, deterministic"

        local wall=0
        local entry
        while read -r entry _
        do
            wall=$((wall + entry))
        done < "${target_dir}/wall.txt"
        echo "  wall time: ${wall} s across $((seeded_count + seedless_count)) run(s), per run in ${target_dir}/wall.txt"
        echo "             (capacity planning only -- NOT a cost measure and not comparable across runs)"

        if [ -n "${waiver}" ] && [ "${findings}" -gt 0 ]; then
            echo "  verdict:   KNOWN RED, ${findings} finding(s) -- ${waiver}"
            echo "KNOWN_RED" > "${campaign_dir}/${name}.status"
        elif [ -n "${waiver}" ]; then
            echo "  verdict:   STALE WAIVER, no finding across the whole campaign."
            echo "             Listed known-red (${waiver}) and did not reproduce at this"
            echo "             execution count. Re-triage at a larger one before removing the"
            echo "             line from tests/fuzz/known_red.txt; this is not a failure."
            echo "STALE_WAIVER" > "${campaign_dir}/${name}.status"
        elif [ "${findings}" -gt 0 ]; then
            echo "  verdict:   FINDINGS, ${findings}"
            echo "FINDINGS" > "${campaign_dir}/${name}.status"
        else
            echo "  verdict:   CLEAN"
            echo "CLEAN" > "${campaign_dir}/${name}.status"
        fi
        echo
    } > "${report}" 2>&1

    return 0
}

target_names=()
if [ "${#requested_targets[@]}" -gt 0 ]; then
    for name in "${requested_targets[@]}"
    do
        if [ ! -x "${fuzz_bin_dir}/${name}" ]; then
            echo "ERROR: no target binary '${name}' under ${fuzz_bin_dir}." >&2
            exit 1
        fi
        target_names+=("${name}")
    done
else
    for fuzz in "${fuzz_bin_dir}"/fuzz_*
    do
        [ -x "${fuzz}" ] || continue
        [ -f "${fuzz}" ] || continue
        target_names+=("$(basename "${fuzz}")")
    done
fi

if [ "${#target_names[@]}" -eq 0 ]; then
    echo "ERROR: no fuzz target binaries found under ${fuzz_bin_dir}." >&2
    exit 1
fi

seed_count=0
for seed in ${seeds}
do
    seed_count=$((seed_count + 1))
done

echo "Build directory:  ${build_abs}"
echo "Campaign output:  ${campaign_dir}"
echo "Targets:          ${#target_names[@]}"
echo "Per run:          -runs=${runs} executions, -max_len=${max_len}"
echo "Seeded runs:      ${seed_count} per target, seeds ${seeds}"
echo "Seedless runs:    ${seedless_runs} per target"
echo "Concurrency:      ${jobs} target(s) at a time"
echo "Waiver list:      ${known_red_file}"
echo

campaign_started="$(date +%s%N)"

running=0
for name in "${target_names[@]}"
do
    sweep_target "${name}" &
    running=$((running + 1))
    if [ "${running}" -ge "${jobs}" ]; then
        # A worker that dies must not take the campaign down with it: hours of
        # measurement on the other targets would go with it, and the loss would
        # be silent. It is caught by the status check below instead, which
        # reports the target by name rather than leaving a gap in the report.
        wait -n || true
        running=$((running - 1))
    fi
done
wait || true

campaign_finished="$(date +%s%N)"

incomplete=()
for name in "${target_names[@]}"
do
    [ -f "${campaign_dir}/${name}.status" ] || incomplete+=("${name}")
done
if [ "${#incomplete[@]}" -gt 0 ]; then
    echo "ERROR: the sweep did not complete for: ${incomplete[*]}" >&2
    echo "  Their output is under ${campaign_dir}. No posture is established for them." >&2
    exit 1
fi

clean=0
findings_targets=()
known_red_targets=()
stale_targets=()

for name in "${target_names[@]}"
do
    cat "${campaign_dir}/${name}.report"
    case "$(cat "${campaign_dir}/${name}.status")" in
        CLEAN)        clean=$((clean + 1)) ;;
        FINDINGS)     findings_targets+=("${name}") ;;
        KNOWN_RED)    known_red_targets+=("${name}") ;;
        STALE_WAIVER) stale_targets+=("${name}") ;;
    esac
done

echo "--- campaign summary ---"
echo "Targets swept:     ${#target_names[@]}"
echo "Runs per target:   ${seed_count} seeded + ${seedless_runs} seedless"
echo "Executions:        $(( (seed_count + seedless_runs) * runs )) per target, \
$(( (seed_count + seedless_runs) * runs * ${#target_names[@]} )) across the campaign"
echo "Campaign wall:     $(elapsed_seconds "${campaign_started}" "${campaign_finished}") s \
(capacity planning only -- NOT a cost measure)"
echo "Clean:             ${clean}"

if [ "${#known_red_targets[@]}" -gt 0 ]; then
    echo "Known red:         ${known_red_targets[*]}"
fi
if [ "${#stale_targets[@]}" -gt 0 ]; then
    echo "Stale waivers:     ${stale_targets[*]} -- re-triage before deleting the line"
fi
if [ "${#findings_targets[@]}" -gt 0 ]; then
    echo "Findings:          ${findings_targets[*]}"
fi

echo "Artifacts kept under ${campaign_dir}; every finding above names the input that reproduces it."

# The suite's state is never stated without naming what is excluded from it. A
# waived target is excluded from the clean count by construction, so a campaign
# that reports every other target clean has established nothing about the suite
# until the exclusion is named beside the count.
if [ "${#known_red_targets[@]}" -gt 0 ]; then
    echo "Clean across ${clean} of ${#target_names[@]} target(s), with ${known_red_targets[*]} \
excluded by waiver and still producing findings."
elif [ "${#stale_targets[@]}" -gt 0 ]; then
    echo "Clean across ${clean} of ${#target_names[@]} target(s), with ${stale_targets[*]} \
excluded by waiver and producing none at this execution count."
fi

if [ "${#findings_targets[@]}" -gt 0 ]; then
    exit 1
fi
