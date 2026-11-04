#ifndef HPP_GUARD_CTRLPP_MPC_NLOPT_TYPES_H
#define HPP_GUARD_CTRLPP_MPC_NLOPT_TYPES_H

#include <cstdint>

namespace ctrlpp
{

// Raw mma and ccsaq reject equality constraints. The auglag variants run them inside NLOPT_AUGLAG_EQ, whose outer
// penalty absorbs the equalities while the inequalities pass through to the inner solver.
enum class nlopt_algorithm : std::uint8_t
{
    slsqp,
    mma,
    cobyla,
    isres,
    auglag_mma,
    auglag_ccsaq,
    ccsaq
};

template <typename Scalar>
struct nlopt_settings
{
    nlopt_algorithm algorithm{nlopt_algorithm::slsqp};
    Scalar ftol_rel{Scalar{1e-6}};
    Scalar xtol_rel{Scalar{1e-6}};
    std::int32_t max_eval{500};
    Scalar max_time{Scalar{0}};
    Scalar constraint_tol{Scalar{1e-8}};
};

}

#endif
