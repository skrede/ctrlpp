#ifndef HPP_GUARD_CTRLPP_CONTROL_DARE_TYPES_H
#define HPP_GUARD_CTRLPP_CONTROL_DARE_TYPES_H

#include "ctrlpp/util/concepts.h"

#include <Eigen/Core>

#include <cstddef>

namespace ctrlpp
{

// Branch on whether a value is present and read the enumerator as a diagnostic: below the precision's
// resolution boundary which refusal a problem gets depends on the build, not on the input.
// @cite laub1979 -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979, Sec. III
enum class dare_error
{
    non_stabilizable, // fewer than n symplectic eigenvalues inside the unit disk
    non_finite_input, // A, B, Q or R holds NaN or Inf
    singular_a,       // A is rank-deficient, so the pencil's A^{-T} does not exist
    singular_r,       // the equilibrated R is rank-deficient, so the pencil's R^{-1} does not exist
    // U11 is numerically singular: either no stabilizing solution exists or the invariant subspace is
    // too ill-conditioned to separate, and the two are not distinguished.
    singular_u11,
    non_psd_solution, // the extracted P is not positive semi-definite to rounding
    schur_failed,     // the real Schur factorization did not converge
    // A finite, admissible problem whose answer could not be verified at the scalar type's precision.
    arithmetic_limit,
};

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU>
struct dare_result
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

    Eigen::Matrix<Scalar, int(NX), int(NX)> P;
    // The gain the acceptance decision was made on, formed at the equilibrated scale; re-forming it from
    // P is not equivalent. The cross-weight overload returns the caller's gain, K + R^{-1} N'.
    Eigen::Matrix<Scalar, int(NU), int(NX)> K;
    Scalar subspace_separation; // LAPACK DTRSEN SEP analogue; near 1 is well-conditioned
    bool reorder_complete;      // false is LAPACK DTRSEN INFO = 1: a diagnostic, not an error
};

}

#endif
