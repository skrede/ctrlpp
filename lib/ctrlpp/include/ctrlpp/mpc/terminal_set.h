#ifndef HPP_GUARD_CTRLPP_MPC_TERMINAL_SET_H
#define HPP_GUARD_CTRLPP_MPC_TERMINAL_SET_H

/// @brief Terminal constraint set types (ellipsoidal and polytopic) for MPC stability.
///
/// @cite mayne2000 -- Mayne et al., "Constrained model predictive control: Stability and optimality", 2000
/// @cite rawlings2017 -- Rawlings et al., "Model Predictive Control: Theory, Computation, and Design", 2017

#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <cstddef>
#include <variant>

namespace ctrlpp {

// Ellipsoidal set: {x : x'Px <= alpha}
// P must be positive definite, alpha > 0.
template<typename Scalar, std::size_t NX>
struct ellipsoidal_set {
    Matrix<Scalar, NX, NX> P;
    Scalar alpha;
};

// Polytopic set in halfplane representation: {x : Hx <= h}
// Dynamic row count (halfplane count varies), fixed NX columns.
template<typename Scalar, std::size_t NX>
struct polytopic_set {
    Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)> H;
    Eigen::VectorX<Scalar> h;
};

// Terminal set: either ellipsoidal or polytopic.
template<typename Scalar, std::size_t NX>
using terminal_set = std::variant<ellipsoidal_set<Scalar, NX>, polytopic_set<Scalar, NX>>;

}

#endif
