#ifndef HPP_GUARD_CTRLPP_DETAIL_NUMERICAL_DIFF_H
#define HPP_GUARD_CTRLPP_DETAIL_NUMERICAL_DIFF_H

// Central differences with the step h = cbrt(eps) * max(1, |x_j|): minimizing truncation O(h^2) plus
// round-off O(eps / h) over h gives h ~ eps^(1/3), with total error ~ eps^(2/3). Each stencil rounds
// h to the nearest representable difference (temp = x_j + h; h = temp - x_j), so the applied
// perturbation is exact in floating point.
// @cite fornberg1988 -- Fornberg, "Generation of finite difference formulas on arbitrarily spaced grids", Math. Comp. 51:184, 1988
// @cite nocedal2006  -- Nocedal & Wright, "Numerical Optimization", 2nd ed., 2006, Sec. 8.1 (central-scheme FD step h = eps^(1/3)*scale; representable-h rounding p. 196)

#include "ctrlpp/detail/numerical_jacobians.h"

#include <span>
#include <cmath>
#include <limits>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <functional>

namespace ctrlpp::detail
{

// `scratch` holds at least z.size() elements; real-time callers keep it for the problem's lifetime so
// no evaluation allocates.
template <typename Scalar>
void finite_diff_gradient(const std::function<Scalar(std::span<const Scalar>)>& f, std::span<const Scalar> z, std::span<Scalar> grad, std::span<Scalar> scratch)
{
    const auto step_scale = std::cbrt(std::numeric_limits<Scalar>::epsilon());
    const auto n = z.size();

    std::copy(z.begin(), z.end(), scratch.begin());

    for(std::size_t j = 0; j < n; ++j)
    {
        const Scalar h_raw = step_scale * std::max(Scalar{1}, std::abs(z[j]));
        const Scalar temp = z[j] + h_raw;
        const Scalar h = temp - z[j];
        const Scalar orig = scratch[j];

        scratch[j] = orig + h;
        const Scalar f_plus = f(std::span<const Scalar>{scratch.data(), n});

        scratch[j] = orig - h;
        const Scalar f_minus = f(std::span<const Scalar>{scratch.data(), n});

        grad[j] = (f_plus - f_minus) / (Scalar{2} * h);
        scratch[j] = orig;
    }
}

// Allocates its scratch on every call; real-time paths use the overload above.
template <typename Scalar>
void finite_diff_gradient(const std::function<Scalar(std::span<const Scalar>)>& f, std::span<const Scalar> z, std::span<Scalar> grad)
{
    std::vector<Scalar> scratch(z.size());
    finite_diff_gradient<Scalar>(f, z, grad, std::span<Scalar>{scratch.data(), scratch.size()});
}

template <typename Scalar>
void evaluate_central_pair(
        const std::function<void(std::span<const Scalar>, std::span<Scalar>)>& c, std::span<Scalar> z, std::size_t j, Scalar h, std::span<Scalar> c_plus, std::span<Scalar> c_minus)
{
    const Scalar orig = z[j];

    z[j] = orig + h;
    c(std::span<const Scalar>{z}, c_plus);

    z[j] = orig - h;
    c(std::span<const Scalar>{z}, c_minus);

    z[j] = orig;
}

template <typename Scalar>
void finite_diff_jacobian(const std::function<void(std::span<const Scalar>, std::span<Scalar>)>& c, std::int32_t n_constraints, std::span<const Scalar> z, std::span<Scalar> jac)
{
    const auto step_scale = std::cbrt(std::numeric_limits<Scalar>::epsilon());
    const auto n = z.size();
    const auto m = static_cast<std::size_t>(n_constraints);

    std::vector<Scalar> z_mut(z.begin(), z.end());
    std::vector<Scalar> c_plus(m);
    std::vector<Scalar> c_minus(m);

    for(std::size_t j = 0; j < n; ++j)
    {
        const Scalar h_raw = step_scale * std::max(Scalar{1}, std::abs(z[j]));
        const Scalar temp = z[j] + h_raw;
        const Scalar h = temp - z[j];

        evaluate_central_pair<Scalar>(c, z_mut, j, h, c_plus, c_minus);
        for(std::size_t i = 0; i < m; ++i)
            jac[i * n + j] = (c_plus[i] - c_minus[i]) / (Scalar{2} * h);
    }
}

}

#endif
