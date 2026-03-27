#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_TRIDIAGONAL_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_TRIDIAGONAL_H

/// @brief Thomas algorithm for tridiagonal and cyclic tridiagonal systems.
///
/// thomas_solve implements the standard Thomas algorithm (forward sweep + back
/// substitution) for tridiagonal systems. cyclic_thomas_solve uses the
/// Sherman-Morrison formula to reduce a cyclic tridiagonal system to two
/// standard tridiagonal solves.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
/// Automatic Machines and Robots", 2009, Sec. 4.4

#include <cassert>
#include <cstddef>
#include <vector>

namespace ctrlpp::detail
{

/// @brief Solve a tridiagonal system using the Thomas algorithm.
///
/// Solves A*x = d in-place where A is tridiagonal with sub-diagonal a,
/// main diagonal b, and super-diagonal c. The solution overwrites d.
/// Both b and d are modified in-place.
///
/// @param a Sub-diagonal coefficients (a[0] unused), size n
/// @param b Main diagonal coefficients (modified in-place), size n
/// @param c Super-diagonal coefficients (c[n-1] unused), size n
/// @param d Right-hand side (overwritten with solution), size n
///
/// @cite biagiotti2009 -- Sec. 4.4, tridiagonal system for cubic spline velocities
template <typename Scalar>
void thomas_solve(std::vector<Scalar> const& a,
                  std::vector<Scalar>& b,
                  std::vector<Scalar> const& c,
                  std::vector<Scalar>& d)
{
    auto const n = b.size();
    assert(n >= 1);
    assert(a.size() == n);
    assert(c.size() == n);
    assert(d.size() == n);

    if (n == 1) {
        d[0] /= b[0];
        return;
    }

    // Forward sweep
    for (std::size_t i = 1; i < n; ++i) {
        auto const w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        d[i] -= w * d[i - 1];
    }

    // Back substitution
    d[n - 1] /= b[n - 1];
    for (std::size_t i = n - 1; i > 0; --i) {
        d[i - 1] = (d[i - 1] - c[i - 1] * d[i]) / b[i - 1];
    }
}

/// @brief Solve a cyclic tridiagonal system using Sherman-Morrison reduction.
///
/// Solves (A + u*v^T)*x = d where the cyclic tridiagonal matrix has corner
/// elements alpha = A[0][n-1] and beta = A[n-1][0]. Uses two standard Thomas
/// solves on copies of the modified system.
///
/// @param a Sub-diagonal coefficients (a[0] unused), size n
/// @param b Main diagonal coefficients (not modified), size n
/// @param c Super-diagonal coefficients (c[n-1] unused), size n
/// @param d Right-hand side (overwritten with solution), size n
/// @param alpha Corner element A[0][n-1]
/// @param beta Corner element A[n-1][0]
///
/// @cite biagiotti2009 -- Sec. 4.4.2, cyclic tridiagonal for periodic BC
template <typename Scalar>
void cyclic_thomas_solve(std::vector<Scalar> const& a,
                         std::vector<Scalar> const& b,
                         std::vector<Scalar> const& c,
                         std::vector<Scalar>& d,
                         Scalar alpha,
                         Scalar beta)
{
    auto const n = b.size();
    assert(n >= 3);

    // gamma = -b[0] (arbitrary nonzero, choosing -b[0] per standard practice)
    auto const gamma = -b[0];

    // Build modified diagonal: b'[0] = b[0] - gamma, b'[n-1] = b[n-1] - alpha*beta/gamma
    std::vector<Scalar> b_mod(b);
    b_mod[0] = b[0] - gamma;
    b_mod[n - 1] = b[n - 1] - alpha * beta / gamma;

    // Build auxiliary vector u: u[0] = gamma, u[1..n-2] = 0, u[n-1] = beta
    // With v = [1, 0, ..., 0, alpha/gamma], u*v^T restores the corner elements:
    //   (0,0): gamma*1 = gamma, (0,n-1): gamma*(alpha/gamma) = alpha
    //   (n-1,0): beta*1 = beta, (n-1,n-1): beta*(alpha/gamma) -- cancelled by b_mod
    std::vector<Scalar> u(n, Scalar{0});
    u[0] = gamma;
    u[n - 1] = beta;

    // Solve A' * y = d (first Thomas solve on copy)
    std::vector<Scalar> b1(b_mod);
    std::vector<Scalar> y(d);
    thomas_solve(a, b1, c, y);

    // Solve A' * z = u (second Thomas solve on copy)
    std::vector<Scalar> b2(b_mod);
    thomas_solve(a, b2, c, u);

    // v = [1, 0, ..., 0, alpha/gamma]
    // v.dot(y) = y[0] + (alpha/gamma) * y[n-1]
    // v.dot(z) = z[0] + (alpha/gamma) * z[n-1]  (z is now u after solve)
    auto const ag = alpha / gamma;
    auto const vy = y[0] + ag * y[n - 1];
    auto const vz = u[0] + ag * u[n - 1];

    // x = y - (v.dot(y) / (1 + v.dot(z))) * z
    auto const factor = vy / (Scalar{1} + vz);
    for (std::size_t i = 0; i < n; ++i) {
        d[i] = y[i] - factor * u[i];
    }
}

} // namespace ctrlpp::detail

#endif
