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
    Scalar max_time{Scalar{0}};
    Scalar constraint_tol{Scalar{1e-8}};
    warm_start_mode warm_start{warm_start_mode::curvature};
};

struct argmin_slsqp
{
    using algorithm = argmin::kraft_slsqp_policy<>;
};

struct argmin_nw_sqp
{
    using algorithm = argmin::nw_sqp_policy<>;
};

// Compile-time-N form of argmin_nw_sqp for the allocation-free static path: it
// resolves argmin::nw_sqp_policy<NV>, whose state_type
// sizes its decision-vector buffers (x, g, bounds, the QP working set) with
// fixed-size Eigen types when NV is a positive compile-time bound. The default
// NV == Eigen::Dynamic yields nw_sqp_policy<> — byte-identical to argmin_nw_sqp.
// The constraint bound M stays derived internally by argmin (dynamic here).
template <int NV = Eigen::Dynamic>
struct argmin_nw_sqp_static
{
    using algorithm = argmin::nw_sqp_policy<NV>;
};

// Rebind a ctrlpp policy binding's argmin algorithm to a compile-time decision
// dimension N when the underlying argmin policy exposes a rebind<N> hook (the
// SQP families do). Policies without the hook are returned unchanged — they have
// no compile-time-N form and are only ever used on the dynamic path. argmin_solver
// uses this to select argmin's compile-time-N step_budget_solver on the static
// path while leaving the dynamic default (N == Eigen::Dynamic) byte-identical:
// for nw_sqp, rebind<Eigen::Dynamic> is nw_sqp_policy<Eigen::Dynamic>, the same
// type as the default-argument nw_sqp_policy<>.
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
    argmin::gradient_tolerance_criterion,
    argmin::objective_tolerance_rel_criterion,
    argmin::step_tolerance_rel_criterion,
    argmin::stall_tolerance_criterion>;

}

#endif
