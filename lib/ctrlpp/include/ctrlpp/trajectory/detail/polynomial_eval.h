#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_POLYNOMIAL_EVAL_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_POLYNOMIAL_EVAL_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
// Machines and Robots", 2009, Sec. 2.1

#include "ctrlpp/types.h"

#include <array>
#include <cstddef>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t ND, std::size_t N>
Vector<Scalar, ND> horner_eval(const std::array<Vector<Scalar, ND>, N>& c, Scalar tau)
{
    static_assert(N >= 1, "polynomial must have at least one coefficient");
    Vector<Scalar, ND> result = c[N - 1];
    for(std::size_t i = N - 1; i > 0; --i)
    {
        result = (c[i - 1] + tau * result).eval();
    }
    return result;
}

template <typename Scalar, std::size_t ND, std::size_t N>
Vector<Scalar, ND> horner_deriv1(const std::array<Vector<Scalar, ND>, N>& c, Scalar tau)
{
    static_assert(N >= 2, "derivative requires at least 2 coefficients");
    Vector<Scalar, ND> result = static_cast<Scalar>(N - 1) * c[N - 1];
    for(std::size_t i = N - 1; i > 1; --i)
    {
        result = (static_cast<Scalar>(i - 1) * c[i - 1] + tau * result).eval();
    }
    return result;
}

template <typename Scalar, std::size_t ND, std::size_t N>
Vector<Scalar, ND> horner_deriv2(const std::array<Vector<Scalar, ND>, N>& c, Scalar tau)
{
    static_assert(N >= 3, "second derivative requires at least 3 coefficients");
    Vector<Scalar, ND> result = static_cast<Scalar>((N - 1) * (N - 2)) * c[N - 1];
    for(std::size_t i = N - 1; i > 2; --i)
    {
        result = (static_cast<Scalar>((i - 1) * (i - 2)) * c[i - 1] + tau * result).eval();
    }
    return result;
}

template <typename Scalar, std::size_t ND, std::size_t N>
Vector<Scalar, ND> horner_deriv3(const std::array<Vector<Scalar, ND>, N>& c, Scalar tau)
{
    static_assert(N >= 4, "third derivative requires at least 4 coefficients");
    Vector<Scalar, ND> result = static_cast<Scalar>((N - 1) * (N - 2) * (N - 3)) * c[N - 1];
    for(std::size_t i = N - 1; i > 3; --i)
    {
        result = (static_cast<Scalar>((i - 1) * (i - 2) * (i - 3)) * c[i - 1] + tau * result).eval();
    }
    return result;
}

}

#endif
