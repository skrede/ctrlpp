#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_TRIDIAGONAL_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_TRIDIAGONAL_H

// The Thomas algorithm, and its cyclic variant reduced to two standard solves through the Sherman-Morrison formula.
// Both refuse a system they cannot tell from a singular one with spline_error::unrepresentable_spline.
// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
// Automatic Machines and Robots", 2009, Sec. 4.4

#include "ctrlpp/trajectory/trajectory_types.h"

#include "ctrlpp/expected.h"

#include <cmath>
#include <limits>
#include <vector>
#include <cassert>
#include <cstddef>
#include <cstdint>

namespace ctrlpp::detail
{

// The computed Thomas solution is exact for A + F with |F| <= (4u + O(u^2)) |L^||U^| (@cite higham1990, eqs. (2.5)-(2.7)), so a smaller pivot cannot be told from zero.
// With b_k the reduced pivot and l_k = a_k / b_(k-1), (|L^||U^|)_kk = |l_k c_(k-1)| + |b_k|; on rows diagonally dominant, as every cubic_spline system is,
// |b_(k-1)| >= |c_(k-1)| by induction, so it is at most |a_k| + |b_k|, which the row scale formed with the reduced pivot bounds.
// The threshold 4 eps = 8u (@cite trefethenbau eq. (13.7)) is twice that first-order bound.
template <typename Scalar>
constexpr Scalar pivot_singular_ulps = Scalar{4};

template <typename Scalar>
Scalar tridiagonal_row_scale(Scalar a_k, Scalar b_k, Scalar c_k)
{
    return std::abs(a_k) + std::abs(b_k) + std::abs(c_k);
}

// Written as a negated comparison so that a NaN value or bound is refused, and strict so that a zero value
// against a zero bound (an all-zero row) is refused too.
template <typename Scalar>
bool indistinguishable_from_zero(Scalar value, Scalar bound)
{
    return !(std::abs(value) > bound);
}

// Checks every pivot once, as it is formed. On refusal b and d are left partially reduced and hold no solution.
template <typename Scalar>
ctrlpp::expected<void, spline_error> thomas_eliminate(std::vector<Scalar> const& a, std::vector<Scalar>& b, std::vector<Scalar> const& c, std::vector<Scalar>& d)
{
    auto const eps = std::numeric_limits<Scalar>::epsilon();
    for(std::size_t i = 0; i < b.size(); ++i)
    {
        if(i > 0)
        {
            auto const w = a[i] / b[i - 1];
            b[i] -= w * c[i - 1];
            d[i] -= w * d[i - 1];
        }
        if(indistinguishable_from_zero(b[i], pivot_singular_ulps<Scalar> * eps * tridiagonal_row_scale(a[i], b[i], c[i])))
        {
            return ctrlpp::unexpected(spline_error::unrepresentable_spline);
        }
    }
    return {};
}

template <typename Scalar>
void thomas_substitute(std::vector<Scalar> const& b, std::vector<Scalar> const& c, std::vector<Scalar>& d)
{
    auto const n = b.size();
    d[n - 1] /= b[n - 1];
    for(std::size_t i = n - 1; i > 0; --i)
    {
        d[i - 1] = (d[i - 1] - c[i - 1] * d[i]) / b[i - 1];
    }
}

// Sub-diagonal a (a[0] enters only the row scale), diagonal b, super-diagonal c (c[n-1] likewise); the solution overwrites d.
// b is overwritten with the reduced pivots, and neither b nor d holds a solution on refusal.
// @cite biagiotti2009 -- Sec. 4.4, tridiagonal system for cubic spline velocities
template <typename Scalar>
ctrlpp::expected<void, spline_error> thomas_solve(std::vector<Scalar> const& a, std::vector<Scalar>& b, std::vector<Scalar> const& c, std::vector<Scalar>& d)
{
    assert(!b.empty());
    assert(a.size() == b.size());
    assert(c.size() == b.size());
    assert(d.size() == b.size());

    auto eliminated = thomas_eliminate(a, b, c, d);
    if(eliminated)
    {
        thomas_substitute(b, c, d);
    }
    return eliminated;
}

// The diagonal of A' = A - u*v^T for u = [gamma, 0, ..., 0, beta] and v = [1, 0, ..., 0, alpha/gamma], whose outer product carries both corners.
template <typename Scalar>
std::vector<Scalar> sherman_morrison_diagonal(std::vector<Scalar> const& b, Scalar alpha, Scalar beta, Scalar gamma)
{
    std::vector<Scalar> b_mod(b);
    b_mod.front() = b.front() - gamma;
    b_mod.back() = b.back() - alpha * beta / gamma;
    return b_mod;
}

template <typename Scalar>
ctrlpp::expected<void, spline_error> sherman_morrison_solves(
        std::vector<Scalar> const& a, std::vector<Scalar> b_mod, std::vector<Scalar> const& c, std::vector<Scalar>& y, std::vector<Scalar>& z)
{
    std::vector<Scalar> b_z(b_mod);
    if(auto solved = thomas_solve(a, b_mod, c, y); !solved)
    {
        return solved;
    }
    return thomas_solve(a, b_z, c, z);
}

template <typename Scalar>
ctrlpp::expected<void, spline_error> sherman_morrison_correct(std::vector<Scalar> const& y, std::vector<Scalar> const& z, Scalar ag, std::vector<Scalar>& d)
{
    // Four roundings (alpha/gamma, its product with z[n-1], the sum forming v^T z, and 1 + v^T z) per @cite trefethenbau
    // eq. (13.7), against the magnitudes that last sum combines. 1 + v^T z is the factor in det(A' + u v^T) = det(A') (1 + v^T A'^{-1} u),
    // from det(I + x y^T) = 1 + y^T x, so it vanishes exactly when the cyclic matrix is singular.
    constexpr std::int32_t sherman_morrison_rounding_ops = 4;
    auto const weighted_last = ag * z.back();
    auto const denominator = Scalar{1} + (z.front() + weighted_last);
    auto const scale = Scalar{1} + std::abs(z.front()) + std::abs(weighted_last);
    if(indistinguishable_from_zero(denominator, static_cast<Scalar>(sherman_morrison_rounding_ops) * std::numeric_limits<Scalar>::epsilon() * scale))
    {
        return ctrlpp::unexpected(spline_error::unrepresentable_spline);
    }
    auto const factor = (y.front() + ag * y.back()) / denominator;
    for(std::size_t i = 0; i < d.size(); ++i)
    {
        d[i] = y[i] - factor * z[i];
    }
    return {};
}

// The corners are alpha = A[0][n-1] and beta = A[n-1][0]; b is not modified, and d is left unchanged on refusal.
// @cite biagiotti2009 -- Sec. 4.4.2, cyclic tridiagonal for periodic BC
template <typename Scalar>
ctrlpp::expected<void, spline_error> cyclic_thomas_solve(
        std::vector<Scalar> const& a, std::vector<Scalar> const& b, std::vector<Scalar> const& c, std::vector<Scalar>& d, Scalar alpha, Scalar beta)
{
    // The reduction is exact down to n = 2 (a periodic spline through 3 waypoints); at n = 1 the corner is the diagonal.
    assert(b.size() >= 2);

    auto const gamma = -b[0];
    std::vector<Scalar> y(d);
    std::vector<Scalar> z(b.size(), Scalar{0});
    z.front() = gamma;
    z.back() = beta;
    if(auto solved = sherman_morrison_solves(a, sherman_morrison_diagonal(b, alpha, beta, gamma), c, y, z); !solved)
    {
        return solved;
    }
    return sherman_morrison_correct(y, z, alpha / gamma, d);
}

}

#endif
