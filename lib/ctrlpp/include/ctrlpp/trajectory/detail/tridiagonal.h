#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_TRIDIAGONAL_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_TRIDIAGONAL_H

/// @brief Thomas algorithm for tridiagonal and cyclic tridiagonal systems.
///
/// thomas_solve implements the standard Thomas algorithm; cyclic_thomas_solve reduces a cyclic
/// tridiagonal system to two standard solves through the Sherman-Morrison formula. Both refuse a
/// system they cannot tell from a singular one with spline_error::unrepresentable_spline.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
/// Automatic Machines and Robots", 2009, Sec. 4.4

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

// The computed Thomas solution is exact for (A + F) with |F| <= (4u + O(u^2)) |L^||U^|: one rounding per
// entry from the factorization, one from forward substitution, two from back substitution (@cite higham1990,
// eqs. (2.5)-(2.7)), so a pivot below that perturbation cannot be told from zero by anything the solver computes.
// Roundings are counted in units of eps = 2u, the per-operation model of @cite trefethenbau eq. (13.7), so 4 eps
// is twice the first-order bound, conservative toward refusing. The scale is the pivot row's coefficient
// magnitudes: for rows diagonally dominant |L^||U^| <= 3|A| (@cite higham1990, Theorem 4.1), and every system
// cubic_spline forms is strictly diagonally dominant by rows.
template <typename Scalar>
constexpr Scalar pivot_singular_ulps = Scalar{4};

template <typename Scalar>
auto tridiagonal_row_scale(Scalar a_k, Scalar b_k, Scalar c_k) -> Scalar
{
    return std::abs(a_k) + std::abs(b_k) + std::abs(c_k);
}

// Written as a negated comparison so that a NaN value or bound is refused, and strict so that a zero value
// against a zero bound (an all-zero row) is refused too.
template <typename Scalar>
auto indistinguishable_from_zero(Scalar value, Scalar bound) -> bool
{
    return !(std::abs(value) > bound);
}

/// @brief Forward elimination, checking every pivot once as it is formed.
///
/// Overwrites b with the reduced pivots and d with the reduced right-hand
/// side. On refusal both are left partially reduced and hold no solution.
template <typename Scalar>
auto thomas_eliminate(std::vector<Scalar> const& a, std::vector<Scalar>& b, std::vector<Scalar> const& c, std::vector<Scalar>& d)
    -> ctrlpp::expected<void, spline_error>
{
    auto const eps = std::numeric_limits<Scalar>::epsilon();
    for (std::size_t i = 0; i < b.size(); ++i)
    {
        if (i > 0)
        {
            auto const w = a[i] / b[i - 1];
            b[i] -= w * c[i - 1];
            d[i] -= w * d[i - 1];
        }
        if (indistinguishable_from_zero(b[i], pivot_singular_ulps<Scalar> * eps * tridiagonal_row_scale(a[i], b[i], c[i])))
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
    for (std::size_t i = n - 1; i > 0; --i)
    {
        d[i - 1] = (d[i - 1] - c[i - 1] * d[i]) / b[i - 1];
    }
}

/// @brief Solve a tridiagonal system using the Thomas algorithm.
///
/// Solves A*x = d in-place where A is tridiagonal with sub-diagonal a,
/// main diagonal b, and super-diagonal c. The solution overwrites d.
/// Both b and d are modified in-place, and hold no solution on refusal.
///
/// @param a Sub-diagonal coefficients (a[0] unused), size n
/// @param b Main diagonal coefficients (modified in-place), size n
/// @param c Super-diagonal coefficients (c[n-1] unused), size n
/// @param d Right-hand side (overwritten with solution), size n
/// @return spline_error::unrepresentable_spline when a pivot is indistinguishable from zero
///
/// @cite biagiotti2009 -- Sec. 4.4, tridiagonal system for cubic spline velocities
template <typename Scalar>
auto thomas_solve(std::vector<Scalar> const& a, std::vector<Scalar>& b, std::vector<Scalar> const& c, std::vector<Scalar>& d)
    -> ctrlpp::expected<void, spline_error>
{
    assert(!b.empty());
    assert(a.size() == b.size());
    assert(c.size() == b.size());
    assert(d.size() == b.size());

    auto eliminated = thomas_eliminate(a, b, c, d);
    if (eliminated)
    {
        thomas_substitute(b, c, d);
    }
    return eliminated;
}

/// @brief Diagonal of A' = A - u*v^T for u = [gamma, 0, ..., 0, beta] and
/// v = [1, 0, ..., 0, alpha/gamma], whose outer product carries both corners.
template <typename Scalar>
auto sherman_morrison_diagonal(std::vector<Scalar> const& b, Scalar alpha, Scalar beta, Scalar gamma) -> std::vector<Scalar>
{
    std::vector<Scalar> b_mod(b);
    b_mod.front() = b.front() - gamma;
    b_mod.back() = b.back() - alpha * beta / gamma;
    return b_mod;
}

template <typename Scalar>
auto sherman_morrison_solves(std::vector<Scalar> const& a, std::vector<Scalar> b_mod, std::vector<Scalar> const& c, std::vector<Scalar>& y,
                             std::vector<Scalar>& z) -> ctrlpp::expected<void, spline_error>
{
    std::vector<Scalar> b_z(b_mod);
    if (auto solved = thomas_solve(a, b_mod, c, y); !solved)
    {
        return solved;
    }
    return thomas_solve(a, b_z, c, z);
}

template <typename Scalar>
auto sherman_morrison_correct(std::vector<Scalar> const& y, std::vector<Scalar> const& z, Scalar ag, std::vector<Scalar>& d)
    -> ctrlpp::expected<void, spline_error>
{
    // Four roundings (alpha/gamma, its product with z[n-1], the sum forming v^T z, and 1 + v^T z) per @cite trefethenbau
    // eq. (13.7), against the magnitudes that last sum combines. 1 + v^T z is the factor in det(A' + u v^T) = det(A') (1 + v^T A'^{-1} u),
    // from det(I + x y^T) = 1 + y^T x, so it vanishes exactly when the cyclic matrix is singular.
    constexpr std::int32_t sherman_morrison_rounding_ops = 4;
    auto const weighted_last = ag * z.back();
    auto const denominator = Scalar{1} + (z.front() + weighted_last);
    auto const scale = Scalar{1} + std::abs(z.front()) + std::abs(weighted_last);
    if (indistinguishable_from_zero(denominator, static_cast<Scalar>(sherman_morrison_rounding_ops) * std::numeric_limits<Scalar>::epsilon() * scale))
    {
        return ctrlpp::unexpected(spline_error::unrepresentable_spline);
    }
    auto const factor = (y.front() + ag * y.back()) / denominator;
    for (std::size_t i = 0; i < d.size(); ++i)
    {
        d[i] = y[i] - factor * z[i];
    }
    return {};
}

/// @brief Solve a cyclic tridiagonal system using Sherman-Morrison reduction.
///
/// Solves (A' + u*v^T)*x = d where the cyclic tridiagonal matrix has corner
/// elements alpha = A[0][n-1] and beta = A[n-1][0].
///
/// @param a Sub-diagonal coefficients (a[0] unused), size n
/// @param b Main diagonal coefficients (not modified), size n
/// @param c Super-diagonal coefficients (c[n-1] unused), size n
/// @param d Right-hand side (overwritten with solution), size n
/// @param alpha Corner element A[0][n-1]
/// @param beta Corner element A[n-1][0]
/// @return spline_error::unrepresentable_spline when either reduced solve is refused or the
/// Sherman-Morrison denominator is indistinguishable from zero; d is then left unchanged
///
/// @cite biagiotti2009 -- Sec. 4.4.2, cyclic tridiagonal for periodic BC
template <typename Scalar>
auto cyclic_thomas_solve(std::vector<Scalar> const& a, std::vector<Scalar> const& b, std::vector<Scalar> const& c, std::vector<Scalar>& d, Scalar alpha,
                         Scalar beta) -> ctrlpp::expected<void, spline_error>
{
    // The reduction is exact down to n = 2 (a periodic spline through 3 waypoints); at n = 1 the corner is the diagonal.
    assert(b.size() >= 2);

    auto const gamma = -b[0];
    std::vector<Scalar> y(d);
    std::vector<Scalar> z(b.size(), Scalar{0});
    z.front() = gamma;
    z.back() = beta;
    if (auto solved = sherman_morrison_solves(a, sherman_morrison_diagonal(b, alpha, beta, gamma), c, y, z); !solved)
    {
        return solved;
    }
    return sherman_morrison_correct(y, z, alpha / gamma, d);
}

}

#endif
