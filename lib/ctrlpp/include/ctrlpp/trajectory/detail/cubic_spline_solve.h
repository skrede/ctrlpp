#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_CUBIC_SPLINE_SOLVE_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_CUBIC_SPLINE_SOLVE_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
// Automatic Machines and Robots", 2009, Sec. 4.4, eq. (4.10)-(4.11)

#include "ctrlpp/trajectory/detail/tridiagonal.h"
#include "ctrlpp/trajectory/detail/spline_checks.h"

#include "ctrlpp/trajectory/trajectory_types.h"

#include "ctrlpp/expected.h"

#include <array>
#include <cmath>
#include <vector>
#include <cstddef>

namespace ctrlpp::detail
{

// The three solves share one parameter list, unused entries included: with them GCC inlines all three into the constructor,
// and without them it keeps the larger two out of line, which deepens construction's stack.
// The natural end rows follow from M_0 = 0 and M_n = 0.
template <typename Scalar>
inline ctrlpp::expected<void, spline_error> solve_natural(
        std::vector<Scalar> const& h, std::vector<Scalar> const& delta, std::vector<Scalar> const& /*pos*/, std::vector<Scalar>& vel, std::size_t n_pts, std::size_t n)
{
    std::vector<Scalar> a(n_pts, Scalar{0});
    std::vector<Scalar> b(n_pts, Scalar{0});
    std::vector<Scalar> c(n_pts, Scalar{0});
    std::vector<Scalar> d(n_pts, Scalar{0});
    b[0] = Scalar{2} * h[0];
    c[0] = h[0];
    d[0] = Scalar{3} * delta[0] * h[0];
    for(std::size_t i = 1; i < n; ++i)
    {
        a[i] = h[i];
        b[i] = Scalar{2} * (h[i - 1] + h[i]);
        c[i] = h[i - 1];
        d[i] = Scalar{3} * (h[i] * delta[i - 1] + h[i - 1] * delta[i]);
    }
    a[n] = h[n - 1];
    b[n] = Scalar{2} * h[n - 1];
    d[n] = Scalar{3} * delta[n - 1] * h[n - 1];
    auto solved = detail::thomas_solve(a, b, c, d);
    if(!solved)
        return solved;
    for(std::size_t i = 0; i < n_pts; ++i)
        vel[i] = d[i];
    return {};
}

// The known end velocities move to the right-hand side of the first and last interior rows.
template <typename Scalar>
inline void set_clamped_row(std::vector<Scalar>& a, std::vector<Scalar>& b, std::vector<Scalar>& c, std::vector<Scalar>& d, std::vector<Scalar> const& h,
        std::vector<Scalar> const& delta, std::size_t j, std::size_t n, Scalar v0, Scalar vn)
{
    auto const n_int = n - 1;
    auto const i = j + 1;
    a[j] = (j > 0) ? h[i] : Scalar{0};
    b[j] = Scalar{2} * (h[i - 1] + h[i]);
    c[j] = (j < n_int - 1) ? h[i - 1] : Scalar{0};
    d[j] = Scalar{3} * (h[i] * delta[i - 1] + h[i - 1] * delta[i]);
    if(i == 1)
    {
        d[j] -= h[i] * v0;
    }
    if(i == n - 1)
    {
        d[j] -= h[i - 1] * vn;
    }
}

// @cite biagiotti2009 -- Sec. 4.4.1
template <typename Scalar>
inline ctrlpp::expected<void, spline_error> solve_clamped(std::vector<Scalar> const& h, std::vector<Scalar> const& delta, std::vector<Scalar> const& /*pos*/, std::vector<Scalar>& vel,
        std::size_t /*n_pts*/, std::size_t n, Scalar v0, Scalar vn)
{
    vel[0] = v0;
    vel[n] = vn;
    if(n <= 1)
    {
        return {};
    }
    auto const n_int = n - 1;
    std::vector<Scalar> a(n_int, Scalar{0});
    std::vector<Scalar> b(n_int, Scalar{0});
    std::vector<Scalar> c(n_int, Scalar{0});
    std::vector<Scalar> d(n_int, Scalar{0});
    for(std::size_t j = 0; j < n_int; ++j)
    {
        set_clamped_row(a, b, c, d, h, delta, j, n, v0, vn);
    }
    auto solved = detail::thomas_solve(a, b, c, d);
    if(!solved)
        return solved;
    for(std::size_t j = 0; j < n_int; ++j)
        vel[j + 1] = d[j];
    return {};
}

// @cite biagiotti2009 -- Sec. 4.4.2
template <typename Scalar>
inline ctrlpp::expected<void, spline_error> solve_periodic(
        std::vector<Scalar> const& h, std::vector<Scalar> const& delta, std::vector<Scalar> const& /*pos*/, std::vector<Scalar>& vel, std::size_t /*n_pts*/, std::size_t n)
{
    std::vector<Scalar> a(n, Scalar{0});
    std::vector<Scalar> b(n, Scalar{0});
    std::vector<Scalar> c(n, Scalar{0});
    std::vector<Scalar> d(n, Scalar{0});
    b[0] = Scalar{2} * (h[n - 1] + h[0]);
    c[0] = h[n - 1];
    d[0] = Scalar{3} * (h[0] * delta[n - 1] + h[n - 1] * delta[0]);
    for(std::size_t i = 1; i < n - 1; ++i)
    {
        a[i] = h[i];
        b[i] = Scalar{2} * (h[i - 1] + h[i]);
        c[i] = h[i - 1];
        d[i] = Scalar{3} * (h[i] * delta[i - 1] + h[i - 1] * delta[i]);
    }
    a[n - 1] = h[n - 1];
    b[n - 1] = Scalar{2} * (h[n - 2] + h[n - 1]);
    d[n - 1] = Scalar{3} * (h[n - 1] * delta[n - 2] + h[n - 2] * delta[n - 1]);
    auto solved = detail::cyclic_thomas_solve(a, b, c, d, h[0], h[n - 2]);
    if(!solved)
        return solved;
    for(std::size_t i = 0; i < n; ++i)
        vel[i] = d[i];
    vel[n] = vel[0];
    return {};
}

// An underflowed quotient comes back subnormal, short of significand, but both numerators are differences of slopes that
// collinear waypoints cancel to a rounding residual; so a subnormal coefficient is refused only when its term's reach
// across the span (numerator times span) exceeds the resolution of the evaluated position.
template <typename Scalar>
inline bool cubic_term_survived(Scalar numerator, Scalar quotient, Scalar span, Scalar scale)
{
    auto const reach = std::abs(numerator * span);
    return (reach <= spline_position_resolution(scale)) || std::isnormal(quotient);
}

template <typename Scalar>
inline bool cubic_coefficients_representable(std::array<Scalar, 4> const& c, Scalar quadratic_numerator, Scalar cubic_numerator, Scalar span, Scalar scale)
{
    bool const finite = std::isfinite(c[0]) && std::isfinite(c[1]) && std::isfinite(c[2]) && std::isfinite(c[3]);

    return finite && cubic_term_survived(quadratic_numerator, c[2], span, scale) && cubic_term_survived(cubic_numerator, c[3], span, scale);
}

}

#endif
