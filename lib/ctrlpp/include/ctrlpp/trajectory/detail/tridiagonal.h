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

#include <cmath>
#include <limits>
#include <vector>
#include <cassert>
#include <cstddef>

namespace ctrlpp::detail
{

// A tridiagonal pivot is treated as numerically singular when its magnitude
// falls below the local row scale times a small multiple of the machine
// epsilon. Scaling by the row coefficients (not a bare absolute constant) keeps
// the test correct across problem magnitudes and float precisions; the multiple
// is the rounding budget of the elimination step, which forms each reduced pivot
// from a handful of multiply-adds.
// @cite trefethenbau -- Trefethen & Bau, "Numerical Linear Algebra", 1997, Lec. 20-22
template <typename Scalar>
constexpr Scalar pivot_singular_ulps = Scalar{4};

template <typename Scalar>
auto tridiagonal_row_scale(Scalar a_k, Scalar b_k, Scalar c_k) -> Scalar
{
    return std::abs(a_k) + std::abs(b_k) + std::abs(c_k);
}

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

    auto const eps = std::numeric_limits<Scalar>::epsilon();
    auto const singular = [eps](Scalar pivot, Scalar row_scale) {
        return std::abs(pivot) < pivot_singular_ulps<Scalar> * eps * row_scale;
    };

    // Forward sweep with pivot monitoring
    for (std::size_t i = 1; i < n; ++i) {
        if (singular(b[i - 1], tridiagonal_row_scale(a[i - 1], b[i - 1], c[i - 1])))
        {
            // Near-zero pivot: zero out remaining unknowns
            for (std::size_t j = i - 1; j < n; ++j)
                d[j] = Scalar{0};
            return;
        }
        auto const w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        d[i] -= w * d[i - 1];
    }

    // Back substitution
    if (singular(b[n - 1], tridiagonal_row_scale(a[n - 1], b[n - 1], c[n - 1])))
    {
        for (auto& di : d)
            di = Scalar{0};
        return;
    }
    d[n - 1] /= b[n - 1];
    for (std::size_t i = n - 1; i > 0; --i) {
        if (singular(b[i - 1], tridiagonal_row_scale(a[i - 1], b[i - 1], c[i - 1])))
        {
            d[i - 1] = Scalar{0};
            continue;
        }
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
    // The Sherman-Morrison reduction is exact down to n = 2: the rank-1 corner
    // update u*v^T adds alpha and beta onto the super- and sub-diagonal entries,
    // which is precisely the dense periodic matrix when the corners and the
    // off-diagonals coincide. n = 2 arises from a periodic cubic spline with 3
    // waypoints, the smallest periodic configuration. n = 1 degenerates (the
    // corner is the diagonal) and stays excluded.
    assert(n >= 2);

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

}

#endif
