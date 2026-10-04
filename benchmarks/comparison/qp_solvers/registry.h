#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_REGISTRY_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_REGISTRY_H

#include <array>
#include <cstdint>
#include <string_view>

namespace ctrlpp::bench::qp_solvers
{

enum class cohort : std::uint8_t
{
    dense,
    sparse
};

enum class tier : std::uint8_t
{
    mpc,
    high
};

enum class regime : std::uint8_t
{
    fixed,
    closed_loop,
    churning
};

enum class dual_kind : std::uint8_t
{
    native,
    recovered
};

enum class update_api : std::uint8_t
{
    separate,
    fused
};

enum class arm_kind : std::uint8_t
{
    argmin_dense_admm,
    argmin_dense_active_set,
    argmin_sparse_admm,
    osqp,
    qpoases,
    daqp,
    proxqp_dense,
    proxqp_sparse,
    hpipm,
    tinympc
};

/// A report-only "polished" arm is its base solver with polishing forced on at the mpc tier, so it
/// shares the base solver's adapter.
struct arm_spec
{
    cohort group;
    std::string_view id;
    arm_kind kind;
    std::string_view role;
    dual_kind duals;
    bool runs_high;
    bool forced_polish;
    update_api update;
};

inline constexpr std::array<arm_spec, 14> arm_registry{{
    {cohort::dense, "argmin_dense_admm", arm_kind::argmin_dense_admm, "candidate", dual_kind::native, true, false, update_api::fused},
    {cohort::dense, "argmin_dense_active_set", arm_kind::argmin_dense_active_set, "candidate", dual_kind::native, true, false, update_api::fused},
    {cohort::dense, "qpoases", arm_kind::qpoases, "competitor", dual_kind::native, true, false, update_api::fused},
    {cohort::dense, "daqp", arm_kind::daqp, "competitor", dual_kind::native, true, false, update_api::separate},
    {cohort::dense, "proxqp_dense", arm_kind::proxqp_dense, "competitor", dual_kind::native, true, false, update_api::separate},
    {cohort::dense, "osqp", arm_kind::osqp, "competitor", dual_kind::native, true, false, update_api::separate},
    {cohort::dense, "argmin_dense_admm_polished", arm_kind::argmin_dense_admm, "report_only", dual_kind::native, false, true, update_api::fused},
    {cohort::sparse, "argmin_sparse_admm", arm_kind::argmin_sparse_admm, "candidate", dual_kind::native, true, false, update_api::fused},
    {cohort::sparse, "osqp", arm_kind::osqp, "competitor", dual_kind::native, true, false, update_api::separate},
    {cohort::sparse, "proxqp_sparse", arm_kind::proxqp_sparse, "competitor", dual_kind::native, true, false, update_api::separate},
    {cohort::sparse, "hpipm", arm_kind::hpipm, "reference", dual_kind::native, true, false, update_api::separate},
    {cohort::sparse, "tinympc", arm_kind::tinympc, "reference", dual_kind::recovered, true, false, update_api::separate},
    {cohort::sparse, "argmin_sparse_admm_polished", arm_kind::argmin_sparse_admm, "report_only", dual_kind::native, false, true, update_api::fused},
    {cohort::sparse, "osqp_polished", arm_kind::osqp, "report_only", dual_kind::native, false, true, update_api::separate},
}};

inline constexpr std::int32_t repetitions = 5;
inline constexpr std::int32_t steps_per_sequence = 60;
inline constexpr std::int32_t fixed_warmup_resolves = 20;
inline constexpr std::int32_t iteration_cap = 4000;
inline constexpr std::int32_t osqp_adaptive_rho_interval = 25;
inline constexpr double mpc_solver_tolerance = 1e-3;
inline constexpr double high_solver_tolerance = 1e-6;
inline constexpr double reference_solver_tolerance = 1e-10;
inline constexpr double mpc_threshold = 1e-3;
inline constexpr double high_threshold = 1e-6;
inline constexpr double probe_tolerance = 1e-6;
inline constexpr double reference_kkt_max = 1e-9;
inline constexpr double binding_row_tolerance = 1e-6;
inline constexpr double binding_median_fraction_min = 0.10;
inline constexpr double binding_state_row_step_fraction_min = 0.50;

inline constexpr std::array<tier, 2> tiers{tier::mpc, tier::high};
inline constexpr std::array<regime, 3> regimes{regime::fixed, regime::closed_loop, regime::churning};

inline auto threshold(tier level) -> double
{
    return level == tier::mpc ? mpc_threshold : high_threshold;
}

inline auto solver_tolerance(tier level) -> double
{
    return level == tier::mpc ? mpc_solver_tolerance : high_solver_tolerance;
}

inline auto name_of(cohort group) -> std::string_view
{
    return group == cohort::dense ? "dense" : "sparse";
}

inline auto name_of(tier level) -> std::string_view
{
    return level == tier::mpc ? "mpc" : "high";
}

inline auto name_of(regime kind) -> std::string_view
{
    constexpr std::array<std::string_view, 3> names{"fixed", "closed_loop", "churning"};
    return names[static_cast<std::size_t>(kind)];
}

inline auto name_of(update_api api) -> std::string_view
{
    return api == update_api::separate ? "separate" : "fused";
}

inline bool runs_at(const arm_spec& arm, tier level)
{
    return level == tier::mpc || arm.runs_high;
}

}

#endif
