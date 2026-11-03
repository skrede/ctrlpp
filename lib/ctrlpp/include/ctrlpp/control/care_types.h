#ifndef HPP_GUARD_CTRLPP_CONTROL_CARE_TYPES_H
#define HPP_GUARD_CTRLPP_CONTROL_CARE_TYPES_H

#include "ctrlpp/util/concepts.h"

#include <Eigen/Core>

#include <cstddef>

namespace ctrlpp
{

// Below the precision's resolution boundary which refusal an input carries is decided by rounding, so branch on
// has_value(); only non_finite_input and singular_r report a property of the operands at any distance from it.
enum class care_error
{
    // An uncontrollable mode at Re > 0 brings its reflection into the left half-plane, so that pair reaches
    // singular_u11 instead.
    non_lhp_stabilizable, // fewer than n Hamiltonian eigenvalues in the open left half-plane
    non_finite_input,     // A, B, Q, R or the assembled Hamiltonian holds a NaN or infinity
    singular_r,           // R is rank-deficient to the reciprocal-pivot tolerance
    // In exact arithmetic a singular U11 means the pair is not both stabilizable and detectable, but the numerical
    // rank test also catches a well-posed pair whose invariant subspace is badly conditioned.
    // @cite laub1979 -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979, Sec. III
    singular_u11,
    non_psd_solution,        // the extracted P is not positive semidefinite within tolerance
    schur_failed,            // Eigen::RealSchur did not converge (Schur methods only)
    sign_function_stagnated, // the sign-function Newton iteration produced nothing it could verify
    unverified_solution,     // a Schur method's extracted P failed the residual or closed-loop postconditions
};

template <ctrlpp_floating_scalar Scalar, std::size_t NX>
struct care_result
{
    static_assert(NX > 0, "State dimension NX must be positive");

    Eigen::Matrix<Scalar, int(NX), int(NX)> P;
    // The minimum pivot ratio across accepted Schur swaps; quiet NaN for a method with no swap phase.
    Scalar subspace_separation;
    bool reorder_complete;
};

}

#endif
