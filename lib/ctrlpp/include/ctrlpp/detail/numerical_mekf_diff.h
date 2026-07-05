#ifndef HPP_GUARD_CTRLPP_DETAIL_NUMERICAL_MEKF_DIFF_H
#define HPP_GUARD_CTRLPP_DETAIL_NUMERICAL_MEKF_DIFF_H

/// @brief MEKF measurement-Jacobian assembly via central differences on the SO(3) tangent space.
///
/// Rotational columns are formed by perturbing the nominal quaternion on the
/// right by so3::exp(+/- step_scale * e_i); bias columns are perturbed in
/// Euclidean space and rounded to a representable difference. This is the
/// Markley 2003 / Sola 2018 error-state convention used by the multiplicative
/// EKF. The step scale follows the central-difference optimum h ~ eps^(1/3)
/// (see numerical_diff.h).
///
/// @cite markley2003 -- Markley, "Attitude Error Representations for Kalman Filtering", 2003
/// @cite sola2018    -- Sola, Deray & Atchuthan, "A Micro Lie Theory for State Estimation in Robotics", 2018, Sec. 5.1
/// @cite fornberg1988 -- Fornberg, "Generation of finite difference formulas on arbitrarily spaced grids", Math. Comp. 51:184, 1988
/// @cite nocedal2006  -- Nocedal & Wright, "Numerical Optimization", 2nd ed., 2006, Sec. 8.1

#include "ctrlpp/types.h"

#include "ctrlpp/lie/so3.h"

#include <Eigen/Geometry>

#include <cmath>
#include <limits>
#include <cstddef>

namespace ctrlpp::detail
{

/// Compute dh/d(error state) via central differences for h(q, b) -> Vector<NY>.
/// Perturbation for rotational columns: q_perturbed = q * so3::exp(step_scale * e_i)
/// Perturbation for bias columns: b_perturbed[j] += step_scale (rounded to representable)
/// Returns Matrix<Scalar, NY, 3+NB>.
template <typename Scalar, std::size_t NB, std::size_t NY, typename H>
auto numerical_mekf_jacobian(const H& h, const Eigen::Quaternion<Scalar>& q, const Vector<Scalar, NB>& b, Scalar step_scale = std::cbrt(std::numeric_limits<Scalar>::epsilon())) -> Matrix<Scalar, NY, 3 + NB>
{
    constexpr std::size_t NE = 3 + NB;
    Matrix<Scalar, NY, NE> jac;
    const Scalar inv_2step = Scalar{1} / (Scalar{2} * step_scale);

    // Rotational perturbation columns (0..2). The error-state baseline is
    // exactly zero, so the perturbation is already exactly representable.
    for(std::size_t i = 0; i < 3; ++i)
    {
        Vector<Scalar, 3> e_i = Vector<Scalar, 3>::Zero();
        e_i(static_cast<Eigen::Index>(i)) = Scalar{1};

        Vector<Scalar, 3> pe = (step_scale * e_i).eval();
        Vector<Scalar, 3> me = (-step_scale * e_i).eval();
        auto q_plus = (q * so3::exp(pe)).normalized();
        auto q_minus = (q * so3::exp(me)).normalized();

        jac.col(static_cast<Eigen::Index>(i)) = (h(q_plus, b) - h(q_minus, b)) * inv_2step;
    }

    // Bias perturbation columns (3..3+NB-1)
    for(std::size_t j = 0; j < NB; ++j)
    {
        Vector<Scalar, NB> b_plus = b;
        Vector<Scalar, NB> b_minus = b;
        const auto idx = static_cast<Eigen::Index>(j);

        const Scalar temp = b(idx) + step_scale;
        const Scalar delta = temp - b(idx);

        b_plus(idx) += delta;
        b_minus(idx) -= delta;

        jac.col(static_cast<Eigen::Index>(3 + j)) = (h(q, b_plus) - h(q, b_minus)) / (Scalar{2} * delta);
    }

    return jac;
}

}

#endif
