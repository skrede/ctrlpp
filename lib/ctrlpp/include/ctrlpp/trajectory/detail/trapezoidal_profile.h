#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_TRAPEZOIDAL_PROFILE_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_TRAPEZOIDAL_PROFILE_H

#include <cmath>
#include <limits>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar>
inline Scalar positive_frame_sign(Scalar h)
{
    return h >= Scalar{0} ? Scalar{1} : Scalar{-1};
}

// x y / (2 z w) with each operand split into fraction and exponent, so no intermediate product leaves the
// range and only the result rounds into it. z and w are positive and finite.
template <typename Scalar>
inline Scalar half_ratio_of_products(Scalar x, Scalar y, Scalar z, Scalar w)
{
    int x_exponent{};
    int y_exponent{};
    int z_exponent{};
    int w_exponent{};
    Scalar const x_fraction = std::frexp(x, &x_exponent);
    Scalar const y_fraction = std::frexp(y, &y_exponent);
    Scalar const z_fraction = std::frexp(z, &z_exponent);
    Scalar const w_fraction = std::frexp(w, &w_exponent);
    return std::ldexp(Scalar{0.5} * x_fraction * y_fraction / (z_fraction * w_fraction), x_exponent + y_exponent - z_exponent - w_exponent);
}

template <typename Scalar>
inline Scalar half_abs_difference_of_squares(Scalar lhs, Scalar rhs)
{
    if(std::abs(lhs) == std::abs(rhs))
    {
        return Scalar{0};
    }

    Scalar const difference = std::abs(lhs - rhs);
    Scalar const sum = std::abs(lhs + rhs);
    if(!std::isfinite(difference) || !std::isfinite(sum))
    {
        return std::numeric_limits<Scalar>::infinity();
    }
    return half_ratio_of_products(difference, sum, Scalar{1}, Scalar{1});
}

// The power-of-two scaling is exact, where std::hypot's three-operand form promises no scaling.
// @cite blue1978 -- Blue, "A Portable Fortran Program to Find the Euclidean
// Norm of a Vector", 1978; anderson2017 -- Anderson, "Algorithm 978: Safe
// Scaling in the Level 1 BLAS", 2017
template <typename Scalar>
inline Scalar scaled_hypot(Scalar x, Scalar y, Scalar z)
{
    int exponent{};
    std::frexp(std::max({x, y, z}), &exponent);
    Scalar const xs = std::ldexp(x, -exponent);
    Scalar const ys = std::ldexp(y, -exponent);
    Scalar const zs = std::ldexp(z, -exponent);
    return std::ldexp(std::sqrt(xs * xs + ys * ys + zs * zs), exponent);
}

// Outside the normal range the product form of eq. (3.14) reads as satisfied on commands that do not
// satisfy it, so the required acceleration is formed as one scaled quotient instead; velocities above
// half the largest finite value are halved first, exactly.
template <typename Scalar, typename Config>
inline Scalar solve_scaled_trapezoidal_acceleration(Config const& cfg, Scalar abs_h)
{
    if(!(abs_h > Scalar{0}))
    {
        return std::numeric_limits<Scalar>::infinity();
    }
    bool const halve = std::max(std::abs(cfg.v0), std::abs(cfg.v1)) > std::numeric_limits<Scalar>::max() / Scalar{2};
    auto const s = halve ? Scalar{0.5} : Scalar{1};
    auto const required = std::ldexp(half_ratio_of_products(std::abs(s * cfg.v0 - s * cfg.v1), std::abs(s * cfg.v0 + s * cfg.v1), abs_h, Scalar{1}), halve ? 2 : 0);
    if(cfg.a_max >= required)
    {
        return cfg.a_max;
    }
    return std::nextafter(required, std::numeric_limits<Scalar>::infinity());
}

// eq. (3.14) made necessary as well as sufficient: in the positive frame the shape exists exactly when the
// triangular peak reaches the larger boundary velocity, which needs a raise only when v0 + v1 > 0.
// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2008, Sec. 3.2.7, eq. (3.14), p.71
template <typename Scalar, typename Config>
inline Scalar solve_trapezoidal_acceleration(Config const& cfg)
{
    auto const h = cfg.q1 - cfg.q0;
    if(!(positive_frame_sign(h) * (cfg.v0 + cfg.v1) > Scalar{0}))
    {
        return cfg.a_max;
    }
    auto const abs_h = std::abs(h);
    auto const v_diff_sq = half_abs_difference_of_squares(cfg.v0, cfg.v1);
    if(!std::isnormal(v_diff_sq) && std::abs(cfg.v0) != std::abs(cfg.v1))
    {
        return solve_scaled_trapezoidal_acceleration(cfg, abs_h);
    }
    if(cfg.a_max * abs_h >= v_diff_sq)
    {
        return cfg.a_max;
    }
    return std::nextafter(v_diff_sq / abs_h, std::numeric_limits<Scalar>::infinity());
}

// The nominal radicand keeps its own rounding wherever it is a finite normal number; past the top of the
// range it overflows and below the smallest normal it has lost the precision the ramps divide by.
template <typename Scalar>
inline Scalar triangular_peak(Scalar abs_h, Scalar v0, Scalar v1, Scalar a)
{
    auto const radicand = (Scalar{2} * a * abs_h + v0 * v0 + v1 * v1) / Scalar{2};
    if(std::isfinite(radicand) && radicand >= std::numeric_limits<Scalar>::min())
    {
        return std::sqrt(radicand);
    }
    auto const root_two = std::sqrt(Scalar{2});
    return scaled_hypot(std::sqrt(a) * std::sqrt(abs_h), std::abs(v0) / root_two, std::abs(v1) / root_two);
}

// (v_peak - v_from) / a cancels catastrophically when the rounded peak lies close to v_from, so for a
// positive v_from it is rationalized to (h + (v_to^2 - v_from^2) / (2 a)) / (v_peak + v_from).
// @cite goldberg1991 -- Goldberg, "What Every Computer Scientist Should Know
// About Floating-Point Arithmetic", 1991, Sec. 1.4; higham2002 -- Higham,
// "Accuracy and Stability of Numerical Algorithms", 2nd ed., 2002, Sec. 1.8
template <typename Scalar>
inline Scalar triangular_ramp_duration(Scalar abs_h, Scalar v_peak, Scalar v_from, Scalar v_to, Scalar a)
{
    if(!(v_from > Scalar{0}))
    {
        return (v_peak - v_from) / a;
    }
    auto const largest = std::max({v_peak, v_from, std::abs(v_to)});
    auto const s = largest > std::numeric_limits<Scalar>::max() / Scalar{2} ? Scalar{0.5} : Scalar{1};
    auto const peak_sum = s * v_peak + s * v_from;
    auto const cruise_term = abs_h / peak_sum * s;
    auto const swept_term = half_ratio_of_products(s * v_to - s * v_from, s * v_to + s * v_from, a, s * peak_sum);
    auto const duration = cruise_term + swept_term;
    // On the boundary eq. (3.14) raises the acceleration onto, the two terms cancel to a ramp of exactly
    // zero, so a negative residual within their four roundings is that zero.
    constexpr Scalar rounding_ops{4};
    auto const resolution = rounding_ops * std::numeric_limits<Scalar>::epsilon() * cruise_term;
    return duration < Scalar{0} && -duration <= resolution ? Scalar{0} : duration;
}

}

#endif
