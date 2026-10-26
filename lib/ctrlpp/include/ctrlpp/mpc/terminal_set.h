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
#include <cstdint>
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
///  * not_converged           : the construction reached its iteration bound without a
///                              fixed point; no set is returned.
///  * closed_loop_not_stable  : the closed loop A + B K is not strictly stable (the Stein
///                              equation has no positive definite solution), or it or K
///                              is not finite, so no bounded output admissible set exists.
///  * unbounded_terminal_set  : the rows C A_K^k for k < NX leave a direction in which
///                              the set recedes, so the polytope the construction starts
///                              from is unbounded. One-sided bounds that would close the
///                              set only at later steps are refused as well, conservatively.
///  * origin_not_interior     : a finite state or input bound does not strictly straddle
///                              zero, or a bound is NaN, so the origin is not interior.
///  * dimension_unsupported   : the state dimension exceeds the four that vertex
///                              enumeration is restricted to.
///  * certification_failed    : the constructed set failed its own check on its vertices
///                              (bounded, inside Y, mapped into itself by A_K, each within
///                              the construction's rounding bounds), or a vertex enumeration
///                              returned fewer than NX + 1 vertices, too few to describe a
///                              full-dimensional polytope; no set is returned.
enum class terminal_set_error : std::uint8_t
{
    input_zero_not_interior,
    empty_terminal_set,
    dare_failed,
    halfplanes_truncated,
    not_converged,
    closed_loop_not_stable,
    unbounded_terminal_set,
    origin_not_interior,
    dimension_unsupported,
    certification_failed,
};

/// @brief Ellipsoid {x : x^T P x <= alpha}, P positive definite, alpha > 0.
///
/// A QP cannot hold a quadric, so the MPC encodes this set as its inscribed box in P's
/// eigenbasis. That box is not positively invariant under the LQR law and gives no
/// recursive-feasibility guarantee; for an invariant terminal set use the output
/// admissible set (terminal_ingredients_result::invariant_set).
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
