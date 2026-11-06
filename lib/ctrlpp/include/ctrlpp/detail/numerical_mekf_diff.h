#ifndef HPP_GUARD_CTRLPP_DETAIL_NUMERICAL_MEKF_DIFF_H
#define HPP_GUARD_CTRLPP_DETAIL_NUMERICAL_MEKF_DIFF_H

// Central differences on the SO(3) error state: a rotational column perturbs the nominal quaternion
// on the right by so3::exp(+/- step_scale * e_i), the multiplicative-EKF convention; a bias column
// perturbs in Euclidean space by a step rounded to a representable difference.
// @cite markley2003 -- Markley, "Attitude Error Representations for Kalman Filtering", 2003
// @cite sola2018    -- Sola, Deray & Atchuthan, "A micro Lie theory for state estimation in robotics", arXiv:1812.01537v9, 2021, Sec. II-G, Eqs. 41-42
// @cite fornberg1988 -- Fornberg, "Generation of finite difference formulas on arbitrarily spaced grids", Math. Comp. 51:184, 1988
// @cite nocedal2006  -- Nocedal & Wright, "Numerical Optimization", 2nd ed., 2006, Sec. 8.1

#include "ctrlpp/types.h"

#include "ctrlpp/lie/so3.h"

#include <Eigen/Geometry>

#include <cmath>
#include <limits>
#include <cstddef>

namespace ctrlpp::detail
{

// The error-state baseline is exactly zero, so the rotational step needs no rounding.
template <typename Scalar, std::size_t NB, std::size_t NY, typename H>
void mekf_rotation_columns(const H& h, const Eigen::Quaternion<Scalar>& q, const Vector<Scalar, NB>& b, Scalar step_scale, Matrix<Scalar, NY, 3 + NB>& jac)
{
    const Scalar inv_2step = Scalar{1} / (Scalar{2} * step_scale);

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
}

template <typename Scalar, std::size_t NB, std::size_t NY, typename H>
void mekf_bias_columns(const H& h, const Eigen::Quaternion<Scalar>& q, const Vector<Scalar, NB>& b, Scalar step_scale, Matrix<Scalar, NY, 3 + NB>& jac)
{
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
}

template <typename Scalar, std::size_t NB, std::size_t NY, typename H>
Matrix<Scalar, NY, 3 + NB> numerical_mekf_jacobian(
        const H& h, const Eigen::Quaternion<Scalar>& q, const Vector<Scalar, NB>& b, Scalar step_scale = std::cbrt(std::numeric_limits<Scalar>::epsilon()))
{
    Matrix<Scalar, NY, 3 + NB> jac;
    mekf_rotation_columns<Scalar, NB, NY>(h, q, b, step_scale, jac);
    mekf_bias_columns<Scalar, NB, NY>(h, q, b, step_scale, jac);
    return jac;
}

}

#endif
