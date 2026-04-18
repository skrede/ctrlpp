#ifndef HPP_GUARD_CTRLPP_CONTROL_CARE_TYPES_H
#define HPP_GUARD_CTRLPP_CONTROL_CARE_TYPES_H

/// @brief Public types for the continuous-time algebraic Riccati equation solver.
///
/// `care_error` enumerates the structured failure modes a CARE solve can produce;
/// `care_result` carries the solution P plus diagnostic scalars. Together they form
/// the `std::expected<care_result, care_error>` contract of `ctrlpp::care`.

#include <Eigen/Core>

#include <cstddef>
#include <type_traits>

namespace ctrlpp
{

/// @brief Structured failure modes for `care`.
///
///  * non_lhp_stabilisable : fewer than n eigenvalues of the Hamiltonian spectrum
///                           lie in the open left half-plane.
///  * non_finite_input     : A, B, Q, R or the assembled Hamiltonian H contains NaN/Inf.
///  * singular_u11         : the top-left n x n block of the reordered invariant-subspace
///                           basis U is singular; P cannot be extracted.
///  * non_psd_solution     : extracted P is not positive semi-definite within an
///                           epsilon-scaled tolerance.
///  * schur_failed         : `Eigen::RealSchur` did not converge on the Hamiltonian.
enum class care_error
{
    non_lhp_stabilisable,
    non_finite_input,
    singular_u11,
    non_psd_solution,
    schur_failed,
};

/// @brief Solution payload of `care`. Same shape as `dare_result`; see dare_types.h
/// for field semantics.
template <typename Scalar, std::size_t NX>
struct care_result
{
    static_assert(std::is_floating_point_v<Scalar>, "Scalar must be a floating-point type");
    static_assert(NX > 0, "State dimension NX must be positive");

    Eigen::Matrix<Scalar, int(NX), int(NX)> P;
    Scalar                                  subspace_separation;
    bool                                    reorder_complete;
};

}

#endif
