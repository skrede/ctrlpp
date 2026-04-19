#ifndef HPP_GUARD_CTRLPP_MPC_ARGMIN_POLICIES_H
#define HPP_GUARD_CTRLPP_MPC_ARGMIN_POLICIES_H

#include <nablapp/solver/isres_policy.h>
#include <nablapp/solver/bobyqa_policy.h>
#include <nablapp/solver/cobyla_policy.h>
#include <nablapp/solver/lbfgsb_policy.h>
#include <nablapp/solver/nw_sqp_policy.h>
#include <nablapp/solver/byrd_lbfgsb_policy.h>
#include <nablapp/solver/kraft_slsqp_policy.h>
#include <nablapp/solver/filter_slsqp_policy.h>
#include <nablapp/solver/filter_nw_sqp_policy.h>
#include <nablapp/solver/augmented_lagrangian_policy.h>

#include <cstdint>

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
    using algorithm = nablapp::kraft_slsqp_policy<>;
};

struct argmin_nw_sqp
{
    using algorithm = nablapp::nw_sqp_policy<>;
};

struct argmin_filter_slsqp
{
    using algorithm = nablapp::filter_slsqp_policy<>;
};

struct argmin_filter_nw_sqp
{
    using algorithm = nablapp::filter_nw_sqp_policy<>;
};

struct argmin_auglag
{
    using algorithm = nablapp::augmented_lagrangian_policy<>;
};

struct argmin_cobyla
{
    using algorithm = nablapp::cobyla_policy;
};

struct argmin_isres
{
    using algorithm = nablapp::isres_policy<>;
};

struct argmin_lbfgsb
{
    using algorithm = nablapp::lbfgsb_policy<>;
};

struct argmin_byrd_lbfgsb
{
    using algorithm = nablapp::byrd_lbfgsb_policy<>;
};

struct argmin_bobyqa
{
    using algorithm = nablapp::bobyqa_policy<>;
};

}

#endif
