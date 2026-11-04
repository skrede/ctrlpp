#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_DERIVED_TOLERANCE_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_DERIVED_TOLERANCE_H

#include <limits>
#include <cstddef>

namespace ctrlpp
{

// A step of the shared control loop is u = -(K*x)(0), nx multiplies and nx - 1 adds (the sign flip is
// exact), then x = A*x + B*u: nx*nx multiplies and nx*(nx - 1) adds for A*x, nx*nu multiplies and at most
// nx*nu adds for B*u, and nx adds to combine them. Charging nx*nu adds rather than nx*(nu - 1) keeps the
// count free of a subtraction that would wrap at nu = 0.
constexpr double loop_roundings_per_step(std::size_t nx, std::size_t nu)
{
    const double n = static_cast<double>(nx);
    const double p = static_cast<double>(nu);
    return 2.0 * n * n + 2.0 * n * p + 2.0 * n - 1.0;
}

// Assembling the plant and weights rounds each entry of A, B, Q and R once. The lqr_gain -> dare path forms
// G = B*R^-1*B^T and A^-T, assembles the 2*nx-square symplectic matrix, takes its real Schur decomposition
// with the basis accumulated (25*m^3), reorders the diagonal blocks, extracts P = U21*U11^-1, checks P for
// positive semi-definiteness and solves (R + B^T*P*B)*K = B^T*P*A; no matrix on it is wider than
// m = 2*nx + nu. Seven stages (the solves forming G, A^-T, P and K, the symplectic assembly, the check's
// factorization and the products B^T*P*B and B^T*P*A) are each a product (2*m^3), a Householder QR with
// back-substitution (at most 3*m^3) or a Cholesky-type factorization (m^3/3), charged at 4*m^3 each. The
// reorder makes fewer than m*m swaps of 2-by-2 block pairs, each charged at 4*4^3.
// @cite golub2013 -- Golub & Van Loan, "Matrix Computations", 4th ed., 2013, Sec. 7.5 (Francis QR with the basis accumulated, 25*m^3), Sec. 5.2 (Householder QR), Sec. 1.1 (products)
constexpr double riccati_roundings(std::size_t nx, std::size_t nu)
{
    const double n = static_cast<double>(nx);
    const double p = static_cast<double>(nu);
    const double m = static_cast<double>(2 * nx + nu);
    const double assembly = 2.0 * n * n + n * p + p * p;
    const double schur = 25.0 * m * m * m;
    const double stages = 7.0 * 4.0 * m * m * m;
    const double reorder = 4.0 * 64.0 * m * m;
    return assembly + schur + stages + reorder;
}

// Every counted operation is charged one epsilon at the scale, whether or not it rounds, with no
// cancellation credited, so the count bounds the departure rather than describing it. Epsilon is twice the
// unit roundoff u, and k*u/(1 - k*u) <= 2*k*u whenever k*u <= 1/2, so k epsilons cover gamma_k exactly.
// @cite higham2002 -- Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., 2002, Ch. 3 (one unit in the last place per operation, accumulated without cancellation; Lemma 3.1 for gamma_k)
template <class Scalar>
constexpr double counted_departure_bound(double roundings, double scale)
{
    const double eps = static_cast<double>(std::numeric_limits<Scalar>::epsilon());
    return roundings * eps * scale;
}

// A board result and the host's double reference are each within the counted bound of the exact run, so
// the two may part by the sum.
template <class Scalar>
constexpr double paired_departure_bound(double roundings, double scale)
{
    return counted_departure_bound<Scalar>(roundings, scale) + counted_departure_bound<double>(roundings, scale);
}

// The count treats the solve as forward-stable to its operation count and does not carry the Riccati
// equation's condition number, so for an ill-conditioned plant it is not a bound; the cross-validation
// suite, not this function, answers whether the solve is right.
template <class Scalar>
constexpr double gain_departure_bound(std::size_t nx, std::size_t nu, double gain_norm)
{
    return counted_departure_bound<Scalar>(riccati_roundings(nx, nu), gain_norm);
}

}

#endif
