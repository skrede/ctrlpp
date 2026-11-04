#ifndef HPP_GUARD_CTRLPP_MPC_TERMINAL_SET_H
#define HPP_GUARD_CTRLPP_MPC_TERMINAL_SET_H

// @cite mayne2000 -- Mayne et al., "Constrained model predictive control: Stability and optimality", 2000
// @cite rawlings2017 -- Rawlings et al., "Model Predictive Control: Theory, Computation, and Design", 2017
// @cite borrelli2017 -- Borrelli, Bemporad & Morari, "Predictive Control for Linear and Hybrid Systems", 2017, Ch. 12 (terminal sets and stability)

#include "ctrlpp/types.h"

#include "ctrlpp/control/dare_types.h"

#include <Eigen/Dense>

#include <cstddef>
#include <cstdint>
#include <variant>

namespace ctrlpp
{

enum class terminal_set_error : std::uint8_t
{
    input_zero_not_interior,
    empty_terminal_set,
    halfplanes_truncated,
    not_converged,
    closed_loop_not_stable,
    unbounded_terminal_set,
    origin_not_interior,
    dimension_unsupported,
    certification_failed,
    non_finite_input,
    singular_input_weight,
};

using terminal_ingredients_error = std::variant<terminal_set_error, dare_error>;

template <typename Scalar, std::size_t NX>
struct ellipsoidal_set
{
    Matrix<Scalar, NX, NX> P;
    Scalar alpha;
};

template <typename Scalar, std::size_t NX>
struct polytopic_set
{
    Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)> H;
    Eigen::VectorX<Scalar> h;
};

template <typename Scalar, std::size_t NX>
using terminal_set = std::variant<ellipsoidal_set<Scalar, NX>, polytopic_set<Scalar, NX>>;

}

#endif
