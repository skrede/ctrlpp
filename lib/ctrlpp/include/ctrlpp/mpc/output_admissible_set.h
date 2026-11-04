#ifndef HPP_GUARD_CTRLPP_MPC_OUTPUT_ADMISSIBLE_SET_H
#define HPP_GUARD_CTRLPP_MPC_OUTPUT_ADMISSIBLE_SET_H

// @cite gilbert1991 -- Gilbert & Tan, "Linear systems with state and control constraints: the theory and application of maximal output admissible sets", IEEE TAC 36(9), 1991
// @cite ossareh2023 -- Ossareh & Kolmanovsky, "On complexity bounds for the maximal admissible set of linear time-invariant systems", 2023 (a-priori step bound)

#include "ctrlpp/mpc/terminal_set.h"

#include "ctrlpp/mpc/detail/polytope_recession.h"
#include "ctrlpp/mpc/detail/polytope_redundancy.h"
#include "ctrlpp/mpc/detail/output_admissible_rows.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include <Eigen/Dense>

#include <limits>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <algorithm>

namespace ctrlpp
{
namespace detail
{

// The first t whose rows C A_K^t are all implied ends the construction, since O_t = O_{t+1} implies O_t = O_inf; the
// set is returned only once it passes its own vertex check. Every O_t lies inside the bounded O_{NX-1} with the origin
// interior, so an enumeration with fewer than NX + 1 vertices has failed, and over it every implication would hold.
template <typename Scalar, std::size_t NX>
ctrlpp::expected<polytopic_set<Scalar, NX>, terminal_set_error> grow_to_fixed_point(polytopic_set<Scalar, NX> set, const polytopic_set<Scalar, NX>& outputs,
        const Matrix<Scalar, NX, NX>& A_K, Matrix<Scalar, NX, NX> power, const Matrix<Scalar, NX, NX>& X, Scalar factor, std::int64_t max_steps)
{
    auto vertices = enumerate_vertices<Scalar, NX>(set, factor);
    Scalar const cap = std::max(static_cast<Scalar>(NX), admissible_step_bound<Scalar, NX>(X, outputs, vertices));
    for(std::int64_t t = static_cast<std::int64_t>(NX);; ++t)
    {
        if(vertices.size() < NX + 1)
            return ctrlpp::unexpected(terminal_set_error::certification_failed);
        set = minimal_rows<Scalar, NX>(set, vertices, factor);
        if(!append_unimplied_rows<Scalar, NX>(set, outputs, power, vertices, factor))
            break;
        if(static_cast<Scalar>(t) > cap || t > max_steps)
            return ctrlpp::unexpected(terminal_set_error::not_converged);
        power = power * A_K;
        vertices = enumerate_vertices<Scalar, NX>(set, factor);
    }
    if(!certified_admissible<Scalar, NX>(set, outputs, A_K, factor))
        return ctrlpp::unexpected(terminal_set_error::certification_failed);
    return set;
}

}

template <typename Scalar, std::size_t NX, std::size_t NU>
ctrlpp::expected<polytopic_set<Scalar, NX>, terminal_set_error> compute_output_admissible_set(const Matrix<Scalar, NX, NX>& A_K, const Matrix<Scalar, NU, NX>& K,
        const Vector<Scalar, NU>& u_min, const Vector<Scalar, NU>& u_max, const Vector<Scalar, NX>& x_min, const Vector<Scalar, NX>& x_max, Scalar rounding_factor = Scalar{1},
        std::int64_t max_steps = std::numeric_limits<std::int64_t>::max())
{
    static_assert(NX <= 4, "Output admissible set computation restricted to NX <= 4");
    bool const nan_bound = u_min.array().isNaN().any() || u_max.array().isNaN().any() || x_min.array().isNaN().any() || x_max.array().isNaN().any();
    if(!A_K.allFinite() || !K.allFinite() || nan_bound)
        return ctrlpp::unexpected(terminal_set_error::non_finite_input);
    if(!detail::straddles_zero<Scalar, NX>(x_min, x_max) || !detail::straddles_zero<Scalar, NU>(u_min, u_max))
        return ctrlpp::unexpected(terminal_set_error::origin_not_interior);
    auto const X = detail::stein_lyapunov_matrix<Scalar, NX>(A_K);
    if(!X)
        return ctrlpp::unexpected(terminal_set_error::closed_loop_not_stable);
    auto const outputs = detail::finite_outputs<Scalar, NX, NU>(K, u_min, u_max, x_min, x_max);
    polytopic_set<Scalar, NX> set{.H = Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>(0, static_cast<Eigen::Index>(NX)), .h = Eigen::VectorX<Scalar>(0)};
    Matrix<Scalar, NX, NX> power = Matrix<Scalar, NX, NX>::Identity();
    for(std::size_t k = 0; k < NX; ++k, power = power * A_K)
        detail::append_output_rows<Scalar, NX>(set, outputs, power);
    if(set.H.rows() < static_cast<Eigen::Index>(NX) || detail::has_recession_direction<Scalar, NX>(set, rounding_factor))
        return ctrlpp::unexpected(terminal_set_error::unbounded_terminal_set);
    return detail::grow_to_fixed_point<Scalar, NX>(std::move(set), outputs, A_K, power, *X, rounding_factor, max_steps);
}

}

#endif
