#ifndef HPP_GUARD_CTRLPP_MPC_ARGMIN_POLICIES_H
#define HPP_GUARD_CTRLPP_MPC_ARGMIN_POLICIES_H

#include <argmin/solver/isres_policy.h>
#include <argmin/solver/bobyqa_policy.h>
#include <argmin/solver/cobyla_policy.h>
#include <argmin/solver/lbfgsb_policy.h>
#include <argmin/solver/nw_sqp_policy.h>
#include <argmin/solver/byrd_lbfgsb_policy.h>
#include <argmin/solver/kraft_slsqp_policy.h>
#include <argmin/solver/filter_slsqp_policy.h>
#include <argmin/solver/filter_nw_sqp_policy.h>
#include <argmin/solver/augmented_lagrangian_policy.h>

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
    using algorithm = argmin::kraft_slsqp_policy<>;
};

struct argmin_nw_sqp
{
    using algorithm = argmin::nw_sqp_policy<>;
};

struct argmin_filter_slsqp
{
    using algorithm = argmin::filter_slsqp_policy<>;
};

struct argmin_filter_nw_sqp
{
    using algorithm = argmin::filter_nw_sqp_policy<>;
};

struct argmin_auglag
{
    using algorithm = argmin::augmented_lagrangian_policy<>;
};

struct argmin_cobyla
{
    using algorithm = argmin::cobyla_policy;
};

struct argmin_isres
{
    using algorithm = argmin::isres_policy<>;
};

struct argmin_lbfgsb
{
    using algorithm = argmin::lbfgsb_policy<>;
};

struct argmin_byrd_lbfgsb
{
    using algorithm = argmin::byrd_lbfgsb_policy<>;
};

struct argmin_bobyqa
{
    using algorithm = argmin::bobyqa_policy<>;
};

}

#endif
