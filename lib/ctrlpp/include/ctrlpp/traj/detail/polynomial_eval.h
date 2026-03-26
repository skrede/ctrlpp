#ifndef HPP_GUARD_CTRLPP_TRAJ_DETAIL_POLYNOMIAL_EVAL_H
#define HPP_GUARD_CTRLPP_TRAJ_DETAIL_POLYNOMIAL_EVAL_H

/// @brief Horner evaluation helpers for polynomial trajectory computation.
///
/// Evaluates polynomials and their derivatives using Horner's method to avoid
/// explicit power computation and improve numerical stability.
///
/// Given coefficients c[0..N-1], the polynomial is:
///   p(tau) = c[0] + c[1]*tau + c[2]*tau^2 + ... + c[N-1]*tau^(N-1)
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009, Sec. 2.1

#include "ctrlpp/types.h"

#include <array>
#include <cstddef>

namespace ctrlpp::detail
{

/// @brief Evaluate polynomial via Horner's method: c[0] + tau*(c[1] + tau*(...))
template <typename Scalar, std::size_t ND, std::size_t N>
auto horner_eval(const std::array<Vector<Scalar, ND>, N>& c, Scalar tau) -> Vector<Scalar, ND>
{
    static_assert(N >= 1, "polynomial must have at least one coefficient");
    Vector<Scalar, ND> result = c[N - 1];
    for(std::size_t i = N - 1; i > 0; --i)
    {
        result = (c[i - 1] + tau * result).eval();
    }
    return result;
}

/// @brief Evaluate first derivative via Horner: c[1] + tau*(2*c[2] + tau*(3*c[3] + ...))
template <typename Scalar, std::size_t ND, std::size_t N>
auto horner_deriv1(const std::array<Vector<Scalar, ND>, N>& c, Scalar tau) -> Vector<Scalar, ND>
{
    static_assert(N >= 2, "derivative requires at least 2 coefficients");
    Vector<Scalar, ND> result = static_cast<Scalar>(N - 1) * c[N - 1];
    for(std::size_t i = N - 1; i > 1; --i)
    {
        result = (static_cast<Scalar>(i - 1) * c[i - 1] + tau * result).eval();
    }
    return result;
}

/// @brief Evaluate second derivative via Horner: 2*c[2] + tau*(6*c[3] + ...)
template <typename Scalar, std::size_t ND, std::size_t N>
auto horner_deriv2(const std::array<Vector<Scalar, ND>, N>& c, Scalar tau) -> Vector<Scalar, ND>
{
    static_assert(N >= 3, "second derivative requires at least 3 coefficients");
    Vector<Scalar, ND> result = static_cast<Scalar>((N - 1) * (N - 2)) * c[N - 1];
    for(std::size_t i = N - 1; i > 2; --i)
    {
        result = (static_cast<Scalar>((i - 1) * (i - 2)) * c[i - 1] + tau * result).eval();
    }
    return result;
}

/// @brief Evaluate third derivative via Horner: 6*c[3] + tau*(24*c[4] + ...)
template <typename Scalar, std::size_t ND, std::size_t N>
auto horner_deriv3(const std::array<Vector<Scalar, ND>, N>& c, Scalar tau) -> Vector<Scalar, ND>
{
    static_assert(N >= 4, "third derivative requires at least 4 coefficients");
    Vector<Scalar, ND> result = static_cast<Scalar>((N - 1) * (N - 2) * (N - 3)) * c[N - 1];
    for(std::size_t i = N - 1; i > 3; --i)
    {
        result = (static_cast<Scalar>((i - 1) * (i - 2) * (i - 3)) * c[i - 1] + tau * result).eval();
    }
    return result;
}

} // namespace ctrlpp::detail

#endif
