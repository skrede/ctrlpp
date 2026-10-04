#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_DERIVED_TOLERANCE_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_DERIVED_TOLERANCE_H

#include <limits>
#include <cstddef>

namespace ctrlpp {

/// @brief Roundings one step of the shared control loop charges, at nx states
/// and nu inputs.
///
/// A step is `u = -(K*x)(0)` then `x = A*x + B*u`. The first is nx multiplies
/// and nx - 1 adds; the sign flip is exact and is not charged. `A*x` is nx*nx
/// multiplies and nx*(nx - 1) adds. `B*u` is nx*nu multiplies and at most nx*nu
/// adds -- reducing nu columns needs only nx*(nu - 1), and charging nx*nu
/// instead keeps the count free of a subtraction that would wrap at nu = 0.
/// Combining the two products is nx more adds. That sums to
/// 2*nx*nx + 2*nx*nu + 2*nx - 1.
constexpr double loop_roundings_per_step(std::size_t nx, std::size_t nu)
{
    const double n = static_cast<double>(nx);
    const double p = static_cast<double>(nu);
    return 2.0 * n * n + 2.0 * n * p + 2.0 * n - 1.0;
}

/// @brief Roundings the gain design charges, at nx states and nu inputs.
///
/// Assembling the plant and weights rounds each entry of A, B, Q and R once.
/// The design path is `lqr_gain` -> `dare`: form G = B*R^-1*B^T and A^-T,
/// assemble the 2*nx-square symplectic matrix, take its real Schur
/// decomposition with the orthogonal basis accumulated, reorder the diagonal
/// blocks, extract P = U21*U11^-1, check P for positive semi-definiteness, then
/// solve (R + B^T*P*B)*K = B^T*P*A. No matrix on that path has a side longer
/// than 2*nx or nu, so m = 2*nx + nu bounds every stage's dimension.
///
/// The Schur decomposition with the basis accumulated is 25*m^3. Seven other
/// stages -- the solve that forms G, the solve that forms A^-T, the products
/// that assemble the symplectic matrix, the solve that forms P, the
/// factorization that checks it, the products B^T*P*B and B^T*P*A, and the
/// solve that forms K -- are each a matrix product (2*m^3), a Householder QR
/// with its back-substitution (at most 3*m^3) or a Cholesky-type factorization
/// (m^3/3), so each is charged at the largest of those, 4*m^3. The reorder
/// makes fewer than m*m adjacent swaps, each on a fixed pair of 2-by-2 blocks,
/// so each swap is charged at 4*4^3 by that same per-stage count.
///
/// @cite golub2013 -- Golub & Van Loan, "Matrix Computations", 4th ed., 2013, Sec. 7.5 (Francis QR with the basis accumulated, 25*m^3), Sec. 5.2 (Householder QR), Sec. 1.1 (products)
constexpr double riccati_roundings(std::size_t nx, std::size_t nu)
{
    const double n        = static_cast<double>(nx);
    const double p        = static_cast<double>(nu);
    const double m        = static_cast<double>(2 * nx + nu);
    const double assembly = 2.0 * n * n + n * p + p * p;
    const double schur    = 25.0 * m * m * m;
    const double stages   = 7.0 * 4.0 * m * m * m;
    const double reorder  = 4.0 * 64.0 * m * m;
    return assembly + schur + stages + reorder;
}

/// @brief Largest departure a quantity computed at `Scalar` may show from its
/// exact value, when the computation is charged `roundings` operations at
/// `scale`.
///
/// Every counted operation is charged one epsilon at the scale, whether or not
/// it rounds, with no cancellation credited, so the count bounds the departure
/// from above rather than describing it. Epsilon is twice the unit roundoff u,
/// and k*u/(1 - k*u) <= 2*k*u whenever k*u <= 1/2, so charging k epsilons
/// covers the gamma_k accumulation exactly.
///
/// @cite higham2002 -- Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., 2002, Ch. 3 (one unit in the last place per operation, accumulated without cancellation; Lemma 3.1 for gamma_k)
template<class Scalar>
constexpr double counted_departure_bound(double roundings, double scale)
{
    const double eps = static_cast<double>(std::numeric_limits<Scalar>::epsilon());
    return roundings * eps * scale;
}

/// @brief How far a board result at `Scalar` and the host's double reference
/// may part: each is within the counted bound of the exact run, so the two may
/// part by the sum.
template<class Scalar>
constexpr double paired_departure_bound(double roundings, double scale)
{
    return counted_departure_bound<Scalar>(roundings, scale) + counted_departure_bound<double>(roundings, scale);
}

/// @brief Largest departure any entry of a gain designed at `Scalar` may show
/// from the exact gain, for a gain whose entries sum in magnitude to
/// `gain_norm`.
///
/// The count treats the solve as forward-stable to its operation count. It
/// does not carry the Riccati equation's condition number, so for an
/// ill-conditioned plant it is not a bound; the cross-validation suite, not
/// this function, answers whether the solve is right.
template<class Scalar>
constexpr double gain_departure_bound(std::size_t nx, std::size_t nu, double gain_norm)
{
    return counted_departure_bound<Scalar>(riccati_roundings(nx, nu), gain_norm);
}

}

#endif
