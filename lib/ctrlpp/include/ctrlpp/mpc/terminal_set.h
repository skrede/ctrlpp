#ifndef HPP_GUARD_CTRLPP_MPC_TERMINAL_SET_H
#define HPP_GUARD_CTRLPP_MPC_TERMINAL_SET_H

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
