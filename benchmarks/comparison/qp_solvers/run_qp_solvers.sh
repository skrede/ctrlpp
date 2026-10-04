#!/usr/bin/env bash
#
# run_qp_solvers.sh -- build and prove out the cross-solver QP suite, then take the pinned recorded run.
#
# Usage:
#   run_qp_solvers.sh prepare
#   run_qp_solvers.sh measure <destination-directory> [--advisory]
#
# Environment (all optional):
#   QP_SOLVERS_BUILD   build tree (default <repo>/build/qp-solvers)
#   QP_SOLVERS_JOBS    build parallelism (default 4)
#   QP_SOLVERS_CPU     the CPU the measurement is pinned to (default 2; never 0)
#   ARGMIN_CHECKOUT    argmin checkout whose verdict command validates the run, used read-only
#                      (default <repo>/../argmin)
#
# "prepare" is insensitive to processor contention. "measure" contains no build step: a build
# discovered inside the measurement is the failure the split exists to prevent.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"
build_tree="${QP_SOLVERS_BUILD:-${repo_root}/build/qp-solvers}"
build_jobs="${QP_SOLVERS_JOBS:-4}"
pin_cpu="${QP_SOLVERS_CPU:-2}"
argmin_checkout="${ARGMIN_CHECKOUT:-${repo_root}/../argmin}"
driver_dir="${build_tree}/comparison/qp_solvers"

# Spelled out because ninja aborts on this project's dynamic dependency file.
generator="Unix Makefiles"
march="-march=x86-64-v3"

preregistration_commit=595da516fa770279855eff1eb9bbd03d5c70e767
argmin_pin=ca7ec854a02dc90f0b423cec52bda9a4f306c43a
proxsuite_upstream=b93d7778ffc3299d84b5cb0851022a29bf24a596
competitor_pins=(
"osqp 236713ce9a56c182ac3230d52108f952afce1523"
"qpoases 680f18c8ef0018a120e1604b769f056e8368df97"
"daqp 29a8987d8b9a596c1ab6acf84c17d38b068b73bc"
"hpipm 16a8db59f595bda642f460aa0f62e308f4b76db3"
"blasfeo 9628393214623b29e5dd2a8784f7423edbc93e04"
"tinympc 023f36bd5b27267e1a96ed47a774d120b5707f13"
)

repetitions=5
smoke_steps=3
smoke_cells_dense="nx2_nu1_N10,nx4_nu2_N50,dense_mpc_shaped"
smoke_cells_sparse="nx2_nu1_N10,nx4_nu2_N50"
# The paths the benchmark build reads; a change under any of them makes the harness dirty.
dirty_paths=(benchmarks lib cmake CMakeLists.txt CMakePresets.json)
thread_env=(OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1)

die()
{
    printf 'run_qp_solvers: %s\n' "$*" >&2
    exit 1
}

usage()
{
    sed -n '3,17p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

utc_now()
{
    date -u '+%Y-%m-%dT%H:%M:%SZ'
}

configure_command()
{
    printf '%q ' cmake -S "${repo_root}/benchmarks" -B "${build_tree}" -G "${generator}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
        "-DCMAKE_C_FLAGS=${march}" "-DCMAKE_CXX_FLAGS=${march}" \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
        -DCMAKE_DISABLE_FIND_PACKAGE_ct_optcon=ON \
        -DCTRLPP_FETCH_BENCHMARK_DEPS=ON \
        -DCTRLPP_ARGMIN_SOURCE_DIR= \
        "-DCTRLPP_ARGMIN_GIT_TAG=${argmin_pin}" \
        -DCTRLPP_BENCH_QP_SOLVERS=ON
    printf '\n'
}

configure_trees()
{
    mkdir -p "${build_tree}"
    configure_command > "${build_tree}/configure_command.txt"
    eval "$(cat "${build_tree}/configure_command.txt")" 2>&1 | tee "${build_tree}/configure.log"
}

build_tree_targets()
{
    cmake --build "${build_tree}" -j"${build_jobs}" --target qp_solvers_dense qp_solvers_sparse
}

fetched_commit()
{
    git -C "${build_tree}/_deps/$1-src" rev-parse HEAD
}

cache_value()
{
    sed -n "s/^$1:[A-Z]*=//p" "${build_tree}/CMakeCache.txt"
}

verify_pins()
{
    local entry name pin
    [ -z "$(cache_value CTRLPP_ARGMIN_SOURCE_DIR)" ] || die "CTRLPP_ARGMIN_SOURCE_DIR is set in the build cache; a recorded run uses the pinned argmin"
    [ "$(fetched_commit argmin)" = "${argmin_pin}" ] || die "fetched argmin is $(fetched_commit argmin), not ${argmin_pin}"
    for entry in "${competitor_pins[@]}"; do
        read -r name pin <<<"${entry}"
        [ "$(fetched_commit "${name}")" = "${pin}" ] || die "fetched ${name} is $(fetched_commit "${name}"), not ${pin}"
    done
}

verify_targets()
{
    local log="${build_tree}/configure.log"
    grep -q 'BLASFEO effective TARGET=X64_INTEL_HASWELL LA=HIGH_PERFORMANCE MF=PANELMAJ' "${log}" || die "BLASFEO target not confirmed in ${log}"
    grep -q 'HPIPM effective TARGET=AVX' "${log}" || die "HPIPM target not confirmed in ${log}"
    grep -q 'DAQP FASTER_MATH=OFF NATIVE_OPTIMIZATION=OFF' "${log}" || die "DAQP options not confirmed in ${log}"
    grep 'effective TARGET\|DAQP FASTER_MATH' "${log}" | sort -u
}

# One compile command per competitor, from the compilation database of the same build.
effective_flags()
{
    python3 - "${build_tree}/compile_commands.json" <<'PY'
import json, sys
units = {"osqp": "/_deps/osqp-build/", "qpoases": "/_deps/qpoases-build/", "daqp": "/_deps/daqp-build/",
         "blasfeo": "/_deps/blasfeo-build/", "hpipm": "/_deps/hpipm-build/", "tinympc": "/_deps/tinympc-build/",
         "harness (argmin, proxsuite headers)": "/comparison/qp_solvers/CMakeFiles/ctrlpp_qp_solvers.dir/proxqp_arms"}
commands = json.load(open(sys.argv[1]))
for name, marker in units.items():
    entry = next((c for c in commands if marker in c.get("output", c["file"]) or marker in c["command"]), None)
    words = entry["command"].split() if entry else []
    flags = [w for w in words if w.startswith(("-O", "-D", "-m", "-f", "-g", "-std")) and not w.startswith("-DCMAKE")]
    print(f"{name}: {' '.join(flags) if entry else 'no compile command found'}")
PY
}

harness_dirty()
{
    if [ -n "$(git -C "${repo_root}" status --porcelain -- "${dirty_paths[@]}")" ]; then
        printf 'true\n'
    else
        printf 'false\n'
    fi
}

compiler_version()
{
    "$(cache_value CMAKE_C_COMPILER)" --version | head -1
}

release_flags()
{
    printf '%s %s\n' "$(cache_value "CMAKE_$1_FLAGS")" "$(cache_value "CMAKE_$1_FLAGS_RELEASE")"
}

require_binaries()
{
    local cohort
    [ -d "${build_tree}" ] || die "no build tree at ${build_tree}; run 'prepare' first"
    for cohort in dense sparse; do
        [ -x "${driver_dir}/qp_solvers_${cohort}" ] || die "qp_solvers_${cohort} is not built; run 'prepare' first"
    done
}

# Host state as a JSON object, read from sysfs and procfs only; nothing is changed.
host_state()
{
    local affinity="$1" frequency_files="$2"
    python3 - "${pin_cpu}" "${affinity}" "${frequency_files}" <<'PY'
import glob, json, platform, statistics, sys
cpu, affinity, frequency_glob = sys.argv[1], sys.argv[2], sys.argv[3]
def read(path):
    try:
        return open(path).read().strip()
    except OSError:
        return None
base = f"/sys/devices/system/cpu/cpu{cpu}/cpufreq/"
samples = []
for path in sorted(glob.glob(frequency_glob)):
    samples += [int(line.split()[2]) for line in open(path) if len(line.split()) == 3]
model = next((l.split(":", 1)[1].strip() for l in open("/proc/cpuinfo") if l.startswith("model name")), None)
state = {
    "cpu_model": model, "kernel": platform.release(),
    "pstate_driver": read(base + "scaling_driver"), "pstate_mode": read("/sys/devices/system/cpu/amd_pstate/status"),
    "energy_performance_preference": read(base + "energy_performance_preference"),
    "boost_global": read("/sys/devices/system/cpu/cpufreq/boost"),
    "boost_policies": {p.split("/")[-2]: read(p) for p in sorted(glob.glob("/sys/devices/system/cpu/cpufreq/policy*/boost"))},
    "governor": read(base + "scaling_governor"), "scaling_min_freq": read(base + "scaling_min_freq"),
    "scaling_max_freq": read(base + "scaling_max_freq"), "cpuinfo_max_freq": read(base + "cpuinfo_max_freq"),
    "smt_control": read("/sys/devices/system/cpu/smt/control"), "smt_active": read("/sys/devices/system/cpu/smt/active"),
    "pinned_cpu": int(cpu), "affinity": affinity, "perf_event_paranoid": read("/proc/sys/kernel/perf_event_paranoid"),
    "cmdline": read("/proc/cmdline"),
    "scaling_cur_freq": {"median": statistics.median(samples) if samples else None, "min": min(samples, default=None),
                         "max": max(samples, default=None), "samples": len(samples)},
}
print(json.dumps(state, sort_keys=True))
PY
}

# Every condition a fixed clock needs that the host does not meet, one per line.
clock_failures()
{
    local host="$1" counter="$2"
    python3 - "${host}" "${counter}" "${pin_cpu}" <<'PY'
import json, sys
h, counter, cpu = json.loads(sys.argv[1]), sys.argv[2], sys.argv[3]
checks = [
    (h["boost_global"] == "0", f"global boost is {h['boost_global']}, needs 0"),
    (all(v == "0" for v in h["boost_policies"].values()), "a per-policy boost is not 0"),
    (h["governor"] == "performance", f"governor of cpu{cpu} is {h['governor']}, needs performance"),
    (h["scaling_min_freq"] is not None and h["scaling_min_freq"] == h["scaling_max_freq"],
     f"cpu{cpu} frequency bounds {h['scaling_min_freq']}..{h['scaling_max_freq']} are not equal"),
    (h["smt_control"] == "off", f"SMT control is {h['smt_control']}, needs off"),
    (h["affinity"].strip() == cpu and cpu != "0", f"affinity is '{h['affinity']}', needs the single CPU {cpu} (not 0)"),
    (counter == "true", "the instruction counter does not open"),
]
for ok, message in checks:
    if not ok:
        print(message)
PY
}

pinned()
{
    taskset -c "${pin_cpu}" env "${thread_env[@]}" "$@"
}

describe_harness()
{
    local dest="$1"
    pinned "${driver_dir}/qp_solvers_dense" --describe "${dest}/describe.json"
}

counter_opens()
{
    python3 -c 'import json,sys; print("true" if json.load(open(sys.argv[1]))["environment"]["instruction_counter_rdpmc"] else "false")' "$1/describe.json"
}

run_gate()
{
    local dest="$1"
    "$2" "${driver_dir}/qp_solvers_dense" --binding-gate "${dest}/binding_gate.json" | tee "${dest}/binding_gate.txt"
}

refused_families()
{
    python3 -c 'import json,sys; print(",".join(f for f, r in json.load(open(sys.argv[1])).items() if not r["family_passed"]))' "$1/binding_gate.json"
}

run_driver()
{
    local launcher="$1" cohort="$2" rep="$3" dest="$4" run_id="$5" refused="$6"
    shift 6
    local filters=("$@")
    if [ -n "${refused}" ]; then
        filters+=(--refused-families "${refused}")
    fi
    "${launcher}" "${driver_dir}/qp_solvers_${cohort}" --rep "${rep}" --out "${dest}" --run-id "${run_id}" "${filters[@]}" \
        2> "${dest}/${cohort}_rep${rep}.log"
}

compress_records()
{
    local file
    for file in "$1"/*_rep[0-9].csv; do
        gzip -n -9 -f "${file}"
    done
}

fill_json()
{
    local run_id="$1" started="$2" finished="$3" host="$4" clock="$5"
    python3 - "${run_id}" "${started}" "${finished}" "${host}" "${clock}" "$(fetched_commit argmin)" "$(git -C "${repo_root}" rev-parse HEAD)" \
        "$(harness_dirty)" "$(compiler_version)" "$(release_flags C)" "$(release_flags CXX)" "${proxsuite_upstream}" \
        "$(for e in "${competitor_pins[@]}"; do read -r n _ <<<"${e}"; printf '%s=%s ' "${n}" "$(fetched_commit "${n}")"; done)" <<'PY'
import json, sys
(run_id, started, finished, host, clock, argmin, harness, dirty, compiler, c_flags, cxx_flags, proxsuite, pins) = sys.argv[1:]
competitor = dict(item.split("=") for item in pins.split())
competitor["proxsuite"] = proxsuite
print(json.dumps({"run_id": run_id, "started_utc": started, "finished_utc": finished, "host": json.loads(host),
                  "clock_state": clock, "preregistration_commit": "595da516fa770279855eff1eb9bbd03d5c70e767",
                  "argmin_sha": argmin, "competitor_pins": competitor, "ctrlpp_harness_sha": harness,
                  "ctrlpp_harness_dirty": dirty == "true", "compiler_version": compiler,
                  "c_flags": " ".join(c_flags.split()), "cxx_flags": " ".join(cxx_flags.split())}))
PY
}

assemble()
{
    local dest="$1"
    python3 "${script_dir}/assemble_manifest.py" --template "${script_dir}/manifest_template.json" \
        --describe "${dest}/describe.json" --gate "${dest}/binding_gate.json" --run-dir "${dest}" --fill "${dest}/fill.json"
}

verdict()
{
    PYTHONDONTWRITEBYTECODE=1 python3 "${argmin_checkout}/benchmarks/scripts/qp_head_to_head_verdicts.py" "$@"
}

verify_template()
{
    verdict --print-manifest-template | python3 -c '
import json, sys
live, committed = json.load(sys.stdin), json.load(open(sys.argv[1]))
sys.exit(0 if live == committed else "the committed manifest template differs from the pre-registration")' \
        "${script_dir}/manifest_template.json"
}

verify_record()
{
    local dest="$1"
    shift
    verdict --run-dir "${dest}" --expect-preregistration-commit "${preregistration_commit}" "$@"
}

write_environment()
{
    local dest="$1"
    cp "${build_tree}/configure_command.txt" "${dest}/configure_command.txt"
    effective_flags > "${dest}/effective_flags.txt"
    {
        printf 'compiler          %s\n' "$(compiler_version)"
        printf 'c flags           %s\n' "$(release_flags C)"
        printf 'c++ flags         %s\n' "$(release_flags CXX)"
        printf 'harness           %s (dirty %s)\n' "$(git -C "${repo_root}" rev-parse HEAD)" "$(harness_dirty)"
        printf 'argmin            %s\n' "$(fetched_commit argmin)"
        printf 'pinned cpu        %s\n' "${pin_cpu}"
    } | tee "${dest}/environment.txt"
}

run_id_now()
{
    printf '%s-%s\n' "$(date -u '+%Y%m%dT%H%M%SZ')" "$(git -C "${repo_root}" rev-parse HEAD | cut -c1-10)"
}

smoke_pass()
{
    local dest="$1" run_id="$2" refused="$3" cohort cells
    for cohort in dense sparse; do
        cells="${smoke_cells_dense}"
        [ "${cohort}" = "sparse" ] && cells="${smoke_cells_sparse}"
        run_driver env "${cohort}" 0 "${dest}" "${run_id}" "${refused}" --steps "${smoke_steps}" --cell "${cells}" \
            --on-probe-failure continue || die "the ${cohort} smoke driver failed (see ${dest}/${cohort}_rep0.log)"
        grep '^probe ' "${dest}/${cohort}_rep0.log" || true
    done
}

probes_failed()
{
    grep -h '^probe .* FAILED' "$1"/*_rep0.log || true
}

prepare()
{
    local started_s smoke dest run_id refused started host failed
    started_s=$(date +%s)
    [ -d "${argmin_checkout}/benchmarks/scripts" ] || die "no argmin checkout at ${argmin_checkout} (set ARGMIN_CHECKOUT)"
    verify_template || die "the committed manifest template differs from argmin's --print-manifest-template"
    configure_trees
    build_tree_targets
    verify_pins
    verify_targets
    smoke="$(date +%s)"
    dest="${build_tree}/smoke"
    rm -rf "${dest}"
    mkdir -p "${dest}"
    write_environment "${dest}"
    describe_harness "${dest}"
    run_gate "${dest}" env
    refused="$(refused_families "${dest}")"
    run_id="$(run_id_now)"
    started="$(utc_now)"
    smoke_pass "${dest}" "${run_id}" "${refused}"
    compress_records "${dest}"
    host="$(host_state "$(taskset -pc $$ | sed 's/.*: //')" "${dest}/*_frequency.txt")"
    fill_json "${run_id}" "${started}" "$(utc_now)" "${host}" unfixed > "${dest}/fill.json"
    assemble "${dest}" || die "the manifest could not be assembled"
    printf 'smoke pass: %d s (configure and build: %d s)\n' "$(( $(date +%s) - smoke ))" "$(( smoke - started_s ))"
    failed="$(probes_failed "${dest}")"
    verify_record "${dest}" --schema-only || die "the smoke run ${dest} does not pass the schema check${failed:+ (failed convention probes:
${failed})}"
}

# A run that aborts after its directory exists is still counted in the ledger.
append_ledger()
{
    local ledger="$1" run_id="$2" started="$3" clock="$4" completed="$5"
    printf '{"run_id": "%s", "started_utc": "%s", "finished_utc": "%s", "clock_state": "%s", "completed": %s}\n' \
        "${run_id}" "${started}" "$(utc_now)" "${clock}" "${completed}" >> "${ledger}"
}

refuse_unless_ready()
{
    require_binaries
    [ "$(harness_dirty)" = "false" ] || die "the harness tree is dirty under ${dirty_paths[*]}; measure from a committed tree"
    verify_pins
    [ "${pin_cpu}" != "0" ] || die "the measurement may not be pinned to CPU 0"
    command -v taskset >/dev/null || die "taskset is required to pin the runs"
}

measure_repetitions()
{
    local dest="$1" run_id="$2" refused="$3" rep cohort
    for rep in $(seq 0 $((repetitions - 1))); do
        for cohort in dense sparse; do
            printf '  rep %s  %s\n' "${rep}" "${cohort}"
            run_driver pinned "${cohort}" "${rep}" "${dest}" "${run_id}" "${refused}" || return 1
        done
    done
}

measure()
{
    local dest_root="${1:-}" advisory="${2:-}" scratch affinity host failures clock run_id dest started refused
    [ -n "${dest_root}" ] || die "measure needs a destination directory, and has no default"
    [ -z "${advisory}" ] || [ "${advisory}" = "--advisory" ] || die "unknown argument '${advisory}'"
    refuse_unless_ready
    scratch="$(mktemp -d)"
    describe_harness "${scratch}"
    affinity="$(taskset -c "${pin_cpu}" bash -c 'taskset -pc $$' | sed 's/.*: //')"
    host="$(host_state "${affinity}" "/nonexistent")"
    failures="$(clock_failures "${host}" "$(counter_opens "${scratch}")")"
    rm -rf "${scratch}"
    clock=fixed
    if [ -n "${failures}" ]; then
        clock=unfixed
        printf 'run_qp_solvers: the clock is not fixed:\n%s\n' "${failures}" >&2
        [ -n "${advisory}" ] || die "not measuring; the operator fixes the host state (or passes --advisory for an unpublishable run)"
    fi
    run_id="$(run_id_now)"
    dest="${dest_root}/${run_id}"
    mkdir -p "${dest}"
    started="$(utc_now)"
    trap 'append_ledger "${dest_root}/runs.jsonl" "${run_id}" "${started}" "${clock}" false' ERR
    write_environment "${dest}"
    describe_harness "${dest}"
    run_gate "${dest}" pinned
    refused="$(refused_families "${dest}")"
    measure_repetitions "${dest}" "${run_id}" "${refused}"
    compress_records "${dest}"
    host="$(host_state "${affinity}" "${dest}/*_frequency.txt")"
    fill_json "${run_id}" "${started}" "$(utc_now)" "${host}" "${clock}" > "${dest}/fill.json"
    assemble "${dest}"
    verify_record "${dest}" --schema-only
    trap - ERR
    append_ledger "${dest_root}/runs.jsonl" "${run_id}" "${started}" "${clock}" true
    printf 'measurement complete: %s (clock %s)\n' "${dest}" "${clock}"
}

case "${1:-}" in
    prepare) prepare ;;
    measure) shift; measure "$@" ;;
    -h|--help) usage ;;
    *) usage >&2; exit 2 ;;
esac
