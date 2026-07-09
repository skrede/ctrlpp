#ifndef HPP_GUARD_CTRLPP_MPC_TERMINAL_SET_H
#define HPP_GUARD_CTRLPP_MPC_TERMINAL_SET_H

/// @brief Terminal constraint set types (ellipsoidal and polytopic) for MPC stability.
///
/// @cite mayne2000 -- Mayne et al., "Constrained model predictive control: Stability and optimality", 2000
/// @cite rawlings2017 -- Rawlings et al., "Model Predictive Control: Theory, Computation, and Design", 2017
/// @cite borrelli2017 -- Borrelli, Bemporad & Morari, "Predictive Control for Linear and Hybrid Systems", 2017, Ch. 12 (terminal sets and stability)

#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <cstddef>
#include <variant>

namespace ctrlpp
{

/// @brief Structured failure modes for terminal-set construction.
///
///  * input_zero_not_interior : some input component has u_min(i) >= 0 or
///                              u_max(i) <= 0, so u = 0 (the LQR input at the
///                              origin) is not strictly interior to [u_min, u_max]
///                              and the input-face alpha bound would be meaningless.
///  * empty_terminal_set      : capping alpha against the input and state faces
///                              yields a non-positive or non-finite alpha, i.e. no
///                              consistent bounded ellipsoid exists (the set is empty).
///  * dare_failed             : the terminal cost and gain could not be formed because
///                              the discrete algebraic Riccati solve failed.
///  * halfplanes_truncated    : the polytopic invariant-set halfplane count hit the
///                              resource bound, so the returned representation is a
///                              truncated (and therefore unsound) outer approximation.
///  * not_converged           : backward reachability did not converge within the
///                              iteration budget; the returned set is not invariant.
enum class terminal_set_error
{
    input_zero_not_interior,
    empty_terminal_set,
    dare_failed,
    halfplanes_truncated,
    not_converged,
};

// Ellipsoidal set: {x : x'Px <= alpha}
// P must be positive definite, alpha > 0.
template <typename Scalar, std::size_t NX>
struct ellipsoidal_set
{
    Matrix<Scalar, NX, NX> P;
    Scalar alpha;
};

// Polytopic set in halfplane representation: {x : Hx <= h}
// Dynamic row count (halfplane count varies), fixed NX columns.
template <typename Scalar, std::size_t NX>
struct polytopic_set
{
    Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)> H;
    Eigen::VectorX<Scalar> h;
};

// Terminal set: either ellipsoidal or polytopic.
template <typename Scalar, std::size_t NX>
using terminal_set = std::variant<ellipsoidal_set<Scalar, NX>, polytopic_set<Scalar, NX>>;

}

#endif
