#ifndef HPP_GUARD_CTRLPP_MPC_ARGMIN_POLICIES_H
#define HPP_GUARD_CTRLPP_MPC_ARGMIN_POLICIES_H

#include <argmin/solver/convergence.h>
#include <argmin/solver/lbfgsb_policy.h>
#include <argmin/solver/nw_sqp_policy.h>
#include <argmin/solver/kraft_slsqp_policy.h>
#include <argmin/solver/filter_slsqp_policy.h>
#include <argmin/solver/filter_nw_sqp_policy.h>
#include <argmin/solver/augmented_lagrangian_policy.h>

#include <Eigen/Core>

#include <cstdint>
#include <concepts>
#include <type_traits>

namespace ctrlpp
{

enum class warm_start_mode : std::uint8_t
{
    cold,
    primal_only,
    curvature
};

template <typename Scalar>
struct argmin_settings
{
    Scalar ftol_rel{Scalar{1e-6}};
    Scalar xtol_rel{Scalar{1e-6}};
    int max_eval{500};
    Scalar constraint_tol{Scalar{1e-8}};
    warm_start_mode warm_start{warm_start_mode::curvature};
    /// Threshold on argmin's composite first-order optimality residual -- the
    /// largest of the Lagrangian-stationarity, primal-feasibility,
    /// dual-feasibility and complementarity legs, each in the infinity norm --
    /// below which a solve stops as `nlp_stop_criterion::stationarity`. It is
    /// absolute, not relative to the problem's scale. The default is argmin's
    /// own, so leaving it unset changes nothing.
    ///
    /// @cite byrd1995 -- Byrd, Lu, Nocedal & Zhu, "A Limited Memory Algorithm for Bound Constrained Optimization", SIAM J. Sci. Comput. 16(5):1190-1208, 1995 (the pgtol default of 1e-5)
    Scalar kkt_tol{Scalar{1e-5}};
};

// A max_time of 0 runs the step-budget driver, which reads no clock.
template <typename Scalar>
struct argmin_timed_settings : argmin_settings<Scalar>
{
    Scalar max_time{Scalar{0}};
};

// Only argmin_step_and_time_budget instantiates argmin's time-budgeted driver,
// and with it a clock read; the settings of a step budget have no max_time to set.
struct argmin_step_budget
{
    template <typename Scalar>
    using settings_type = argmin_settings<Scalar>;
};

struct argmin_step_and_time_budget
{
    template <typename Scalar>
    using settings_type = argmin_timed_settings<Scalar>;
};

template <typename Budget>
concept argmin_budget = std::same_as<Budget, argmin_step_budget> || std::same_as<Budget, argmin_step_and_time_budget>;

struct argmin_slsqp
{
    using algorithm = argmin::kraft_slsqp_policy<>;
};

struct argmin_nw_sqp
{
    using algorithm = argmin::nw_sqp_policy<>;
};

// A compile-time NV sizes the policy's decision-vector buffers with fixed-size Eigen types on the allocation-free
// static path; NV == Eigen::Dynamic yields the same type as argmin_nw_sqp's policy.
template <int NV = Eigen::Dynamic>
struct argmin_nw_sqp_static
{
    using algorithm = argmin::nw_sqp_policy<double, NV>;
};

// Policies with a rebind hook are re-dimensioned to the compile-time N and others pass through unchanged; rebinding
// to Eigen::Dynamic leaves the dynamic default's type unchanged.
template <typename Algorithm, int N, typename = void>
struct rebind_argmin_algorithm
{
    using type = Algorithm;
};

template <typename Algorithm, int N>
struct rebind_argmin_algorithm<Algorithm, N,
    std::void_t<typename Algorithm::template rebind<N>>>
{
    using type = typename Algorithm::template rebind<N>;
};

template <typename Algorithm, int N>
using rebind_argmin_algorithm_t = typename rebind_argmin_algorithm<Algorithm, N>::type;

struct argmin_filter_slsqp
{
    using algorithm = argmin::filter_slsqp_policy<>;
};

struct argmin_filter_nw_sqp
{
    using algorithm = argmin::filter_nw_sqp_policy<>;
};

struct argmin_lbfgsb
{
    using algorithm = argmin::lbfgsb_policy<>;
};

template <typename Inner = argmin_lbfgsb>
struct argmin_auglag
{
    using algorithm = argmin::augmented_lagrangian_policy<typename Inner::algorithm>;
};

// Convergence policy used by every ctrlpp argmin_solver instantiation. It
// carries the RELATIVE objective/step criteria (objective_tolerance_rel /
// step_tolerance_rel) so the ctrlpp settings fields named ftol_rel / xtol_rel
// can wire to argmin's `_rel` setters, which are requires-guarded on those
// criteria being present in the convergence tuple. Gradient and stall criteria
// are retained from argmin's default_convergence. This deliberately REPLACES
// default_convergence's ABSOLUTE objective/step criteria with their relative
// counterparts and is a deliberate convergence-behavior change: argmin-backed
// solvers observe re-pinned iterate counts relative to default_convergence.
using argmin_ctrlpp_convergence = argmin::convergence_policy<
    argmin::gradient_tolerance_criterion<>,
    argmin::objective_tolerance_rel_criterion<>,
    argmin::step_tolerance_rel_criterion<>,
    argmin::stall_tolerance_criterion<>>;

}

#endif
