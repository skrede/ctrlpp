#!/usr/bin/env bash
#
# run_window.sh -- build and prove out the benchmark set, then take the pinned
# measurement runs.
#
# Usage:
#   run_window.sh prepare
#   run_window.sh measure <destination-directory>
#
# Everything in "prepare" is insensitive to processor contention and is meant to
# run well before the measurement. "measure" contains no build step at all: a
# build discovered inside the measurement is the failure this split exists to
# prevent, so a missing binary is an error naming the target.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

default_tree="${script_dir}/build-window"
# The nonlinear optimal-control competitor's umbrella assigns to a const member
# at ct/optcon/dms/dms_core/TimeGrid.h:53, which g++ 13 accepts and later
# releases reject, so that one target gets its own tree and its own compiler.
second_tree="${script_dir}/build-window-gcc13"
second_cxx="g++-13"

# Three passes pinned to this one core: both match the runs already archived, so
# a new measurement is comparable with them instead of starting a new baseline.
passes=3
pin_cpu=2

# Spelled out because ninja 1.13.2 aborts on this project's dynamic dependency
# file. Parallelism belongs to the build only; the measurement is serial.
generator="Unix Makefiles"
build_jobs=6

schema_header='"title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total","accuracy_metric","accuracy_value"'

# tree | directory under the tree | target | data rows its shape predicts | window
#
# A dash in the row column marks a benchmark writing its own record format
# rather than the shared twelve-column one; those four also take no smoke
# switch, so preparation runs them whole.
# "hold" marks a target that must keep compiling but whose row is not published,
# so it is given no measurement time.
#
# The two targets whose cost cannot be derived from an earlier archived run are
# listed first, so an unexpectedly long one surfaces while there is still
# measuring time left rather than after the budget is spent.
targets=(
"g comparison/ct bench_nmpc_vs_ct 4 run"
"w comparison/libmpc bench_mpc_vs_libmpc 4 run"
"w internal bench_pid 1 run"
"w internal bench_lqr 4 run"
"w internal bench_dare 108 run"
"w internal bench_kalman 2 run"
"w internal bench_ekf 2 run"
"w internal bench_ukf 2 run"
"w internal bench_mrac 1 run"
"w internal bench_l1 1 run"
"w internal bench_sysid 4 run"
"w internal bench_trajectory 4 run"
"w internal bench_dsp 4 run"
"w internal bench_so3 4 run"
"w internal bench_care_vs_ct_nx16 4 hold"
"w internal bench_care_vs_ct_nx30 4 hold"
"w internal bench_care_methods 81 hold"
"w comparison/ct bench_care_vs_ct 36 run"
"w comparison/ct bench_lqr_continuous_vs_ct 54 run"
"w comparison/ct bench_dare_vs_ct 20 run"
"w comparison/ct bench_lqr_vs_ct 20 hold"
"w comparison/ruckig bench_trajectory_vs_ruckig 10 run"
"w comparison/osqp_eigen bench_qp_vs_osqp_eigen 4 run"
"w comparison/proxsuite bench_qp_vs_proxqp 4 hold"
"w comparison/argmin bench_qp_vs_osqp 24 hold"
"w comparison/argmin bench_jacobian 16 hold"
"w comparison/argmin bench_slsqp 114 hold"
"w comparison/argmin bench_sqp_variants 88 hold"
"w comparison/argmin bench_nmpc 48 hold"
"w comparison/argmin bench_nmhe 8 hold"
"w comparison/argmin bench_schedule 10 hold"
"w comparison/argmin bench_qp_vs_osqp_closedloop - hold"
"w comparison/argmin bench_qp_vs_osqp_churning - run"
"w comparison/argmin bench_qp_periter - hold"
"w comparison/argmin bench_step_budget - hold"
)

die()
{
    printf 'run_window: %s\n' "$*" >&2
    exit 1
}

usage()
{
    cat <<'USAGE'
Usage:
  run_window.sh prepare
      Configure both build trees, build every target the measurement runs, run
      each one under its smoke switch and check the record it writes.

  run_window.sh measure <destination-directory>
      Take the measurement runs into the given directory. Builds nothing. The
      destination is required and has no default: the archive lives outside this
      repository and its path may not appear inside it.
USAGE
}

tree_path()
{
    case "$1" in
        w) printf '%s\n' "${default_tree}" ;;
        g) printf '%s\n' "${second_tree}" ;;
        *) die "unknown build tree '$1'" ;;
    esac
}

# The optimizer's revision is read from the library's own pin rather than
# repeated here, so the measured revision cannot drift from the shipped one.
argmin_pin()
{
    local pin
    pin="$(sed -n 's/^set(CTRLPP_ARGMIN_GIT_TAG "\([^"]*\)".*/\1/p' "${repo_root}/CMakeLists.txt")"
    [ -n "${pin}" ] || die "could not read the optimizer pin from the library build file"
    printf '%s\n' "${pin}"
}

configure_trees()
{
    cmake -S "${script_dir}" -B "${default_tree}" -G "${generator}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCTRLPP_FETCH_BENCHMARK_DEPS=ON \
        -DCTRLPP_ARGMIN_GIT_TAG="$(argmin_pin)" \
        -DCTRLPP_BENCH_CT=ON -DCTRLPP_BENCH_RUCKIG=ON -DCTRLPP_BENCH_LIBMPC=ON \
        -DCTRLPP_BENCH_OSQP_EIGEN=ON -DCTRLPP_BENCH_PROXSUITE=ON \
        -DCTRLPP_BENCH_ARGMIN=ON
    cmake -S "${script_dir}" -B "${second_tree}" -G "${generator}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCTRLPP_FETCH_BENCHMARK_DEPS=ON \
        -DCMAKE_CXX_COMPILER="${second_cxx}" \
        -DCTRLPP_BENCH_CT=ON
}

targets_in_tree()
{
    local wanted="$1" entry tree target
    for entry in "${targets[@]}"; do
        read -r tree _ target _ <<<"${entry}"
        if [ "${tree}" = "${wanted}" ]; then
            printf '%s\n' "${target}"
        fi
    done
}

build_trees()
{
    local tree list
    for tree in w g; do
        mapfile -t list < <(targets_in_tree "${tree}")
        [ "${#list[@]}" -gt 0 ] || continue
        cmake --build "$(tree_path "${tree}")" -j"${build_jobs}" --target "${list[@]}"
    done
}

record_path()
{
    local work="$1" target="$2" rows="$3" candidate
    if [ "${rows}" = "-" ]; then
        candidate="$(find "${work}" -maxdepth 1 -name '*.csv' | sort | head -1)"
    elif [ -f "${work}/${target}_timing.csv" ]; then
        candidate="${work}/${target}_timing.csv"
    else
        candidate="${work}/${target}.csv"
    fi
    [ -n "${candidate}" ] && [ -f "${candidate}" ] || return 1
    printf '%s\n' "${candidate}"
}

verify_record()
{
    local work="$1" target="$2" rows="$3" record observed
    record="$(record_path "${work}" "${target}" "${rows}")" || { printf 'no record written'; return 1; }
    observed=$(( $(wc -l < "${record}") - 1 ))
    if [ "${rows}" = "-" ]; then
        [ "${observed}" -ge 1 ] || { printf 'own record format, empty'; return 1; }
        printf 'own record format, %d data rows' "${observed}"
        return 0
    fi
    [ "$(head -1 "${record}")" = "${schema_header}" ] || { printf 'header is not the shared twelve-column one'; return 1; }
    [ "${observed}" -eq "${rows}" ] || { printf 'wrote %d data rows, its shape predicts %d' "${observed}" "${rows}"; return 1; }
    printf '%d data rows, as its shape predicts' "${observed}"
}

# Each benchmark writes a fixed filename into the working directory, so every
# run gets a directory of its own.
smoke_one()
{
    local tree="$1" rundir="$2" target="$3" rows="$4" work bin
    bin="$(tree_path "${tree}")/${rundir}/${target}"
    work="$(tree_path "${tree}")/smoke/${target}"
    [ -x "${bin}" ] || return 1
    rm -rf "${work}"
    mkdir -p "${work}"
    if [ "${rows}" = "-" ]; then
        ( cd "${work}" && "${bin}" ) >"${work}/stdout.txt" 2>&1
    else
        ( cd "${work}" && "${bin}" --smoke ) >"${work}/stdout.txt" 2>&1
    fi
}

check_one()
{
    local tree="$1" rundir="$2" target="$3" rows="$4" work verdict
    work="$(tree_path "${tree}")/smoke/${target}"
    if ! smoke_one "${tree}" "${rundir}" "${target}" "${rows}"; then
        printf '%-30s FAIL  did not build, or exited nonzero\n' "${target}"
        return 1
    fi
    if ! verdict="$(verify_record "${work}" "${target}" "${rows}")"; then
        printf '%-30s FAIL  %s\n' "${target}" "${verdict}"
        return 1
    fi
    printf '%-30s ok    %s\n' "${target}" "${verdict}"
}

prepare()
{
    configure_trees
    build_trees
    local entry failures=0
    for entry in "${targets[@]}"; do
        # shellcheck disable=SC2086
        check_one ${entry} || failures=$((failures + 1))
    done
    [ "${failures}" -eq 0 ] || die "${failures} of ${#targets[@]} targets are not ready"
    printf '%d targets configured, built, ran, and wrote the record their shape predicts\n' "${#targets[@]}"
}

require_binaries()
{
    local entry tree rundir target window missing=0
    for entry in "${targets[@]}"; do
        read -r tree rundir target _ window <<<"${entry}"
        [ "${window}" = "run" ] || continue
        if [ ! -x "$(tree_path "${tree}")/${rundir}/${target}" ]; then
            printf 'run_window: %s is not built\n' "${target}" >&2
            missing=$((missing + 1))
        fi
    done
    [ "${missing}" -eq 0 ] || die "${missing} target(s) missing; run 'prepare' before the measurement"
}

write_environment()
{
    local dest="$1"
    {
        printf 'taken             %s\n' "$(date -u '+%Y-%m-%d %H:%M:%SZ')"
        printf 'kernel            %s\n' "$(uname -srm)"
        printf 'processor         %s\n' "$(sed -n 's/^model name[[:space:]]*: //p' /proc/cpuinfo | head -1)"
        printf 'governor          %s\n' "$(cat "/sys/devices/system/cpu/cpu${pin_cpu}/cpufreq/scaling_governor" 2>/dev/null || echo unknown)"
        printf 'boost             %s\n' "$(cat /sys/devices/system/cpu/cpufreq/boost 2>/dev/null || echo unknown)"
        printf 'pinning           taskset -c %s\n' "${pin_cpu}"
        printf 'passes            %s\n' "${passes}"
        printf 'default compiler  %s\n' "$(c++ --version | head -1)"
        printf 'second compiler   %s\n' "$("${second_cxx}" --version | head -1)"
        printf 'library revision  %s\n' "$(git -C "${repo_root}" rev-parse --short HEAD)"
        printf 'optimizer pin     %s\n' "$(argmin_pin)"
    } | tee "${dest}/environment.txt"
}

measure_one()
{
    local dest="$1" pass="$2" tree="$3" rundir="$4" target="$5" rows="$6" window="$7" work
    [ "${window}" = "run" ] || return 0
    work="${dest}/pass-${pass}/${target}"
    mkdir -p "${work}"
    printf '  pass %s  %s\n' "${pass}" "${target}"
    ( cd "${work}" && taskset -c "${pin_cpu}" "$(tree_path "${tree}")/${rundir}/${target}" ) \
        >"${work}/stdout.txt" 2>&1 || die "${target} failed in pass ${pass}"
}

measure()
{
    local dest="${1:-}" pass entry
    [ -n "${dest}" ] || die "measure needs a destination directory, and has no default"
    command -v taskset >/dev/null || die "taskset is required to pin the runs"
    mkdir -p "${dest}"
    require_binaries
    write_environment "${dest}"
    for pass in $(seq 1 "${passes}"); do
        for entry in "${targets[@]}"; do
            # shellcheck disable=SC2086
            measure_one "${dest}" "${pass}" ${entry}
        done
    done
    printf 'measurement complete: %d passes under %s\n' "${passes}" "${dest}"
}

case "${1:-}" in
    prepare) prepare ;;
    measure) shift; measure "$@" ;;
    -h|--help) usage ;;
    *) usage >&2; exit 2 ;;
esac
