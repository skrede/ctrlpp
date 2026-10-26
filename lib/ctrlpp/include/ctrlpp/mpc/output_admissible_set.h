#ifndef HPP_GUARD_CTRLPP_MPC_OUTPUT_ADMISSIBLE_SET_H
#define HPP_GUARD_CTRLPP_MPC_OUTPUT_ADMISSIBLE_SET_H

/// @brief Maximal output admissible set of a stable closed loop under state and input bounds.
///
/// @cite gilbert1991 -- Gilbert & Tan, "Linear systems with state and control constraints: the theory and application of maximal output admissible sets", IEEE TAC 36(9), 1991
/// @cite ossareh2023 -- Ossareh & Kolmanovsky, "On complexity bounds for the maximal admissible set of linear time-invariant systems", 2023 (a-priori step bound)

#include "ctrlpp/mpc/terminal_set.h"
#include "ctrlpp/mpc/detail/polytope_recession.h"
#include "ctrlpp/mpc/detail/polytope_redundancy.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <optional>
#include <algorithm>

namespace ctrlpp
{
namespace detail
{

template <typename Scalar, std::size_t N>
bool straddles_zero(const Vector<Scalar, N>& lower, const Vector<Scalar, N>& upper)
{
    return ((lower.array() < Scalar{0}) && (upper.array() > Scalar{0})).all();
}

/// The output map y = [I; K] x with one row per finite bound side, unnormalized; an infinite side constrains nothing.
template <typename Scalar, std::size_t NX, std::size_t NU>
polytopic_set<Scalar, NX> finite_outputs(const Matrix<Scalar, NU, NX>& K, const Vector<Scalar, NU>& u_min, const Vector<Scalar, NU>& u_max,
                                         const Vector<Scalar, NX>& x_min, const Vector<Scalar, NX>& x_max)
{
    polytopic_set<Scalar, NX> outputs{.H = Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>(0, static_cast<Eigen::Index>(NX)), .h = Eigen::VectorX<Scalar>(0)};
    auto add = [&outputs](const polytope_row<Scalar, NX>& row, Scalar lower, Scalar upper)
    {
        if(std::isfinite(upper))
            append_row<Scalar, NX>(outputs, row, upper);
        if(std::isfinite(lower))
            append_row<Scalar, NX>(outputs, -row, -lower);
    };
    for(Eigen::Index j = 0; j < static_cast<Eigen::Index>(NX); ++j)
        add(polytope_row<Scalar, NX>::Unit(j), x_min(j), x_max(j));
    for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NU); ++i)
        add(K.row(i), u_min(i), u_max(i));
    return outputs;
}

/// A zero row is dropped: the origin-interior hypothesis makes its bound positive, so it holds everywhere.
template <typename Scalar, std::size_t NX>
std::optional<std::pair<polytope_row<Scalar, NX>, Scalar>> unit_output_row(const polytopic_set<Scalar, NX>& outputs, Eigen::Index i, const Matrix<Scalar, NX, NX>& power)
{
    polytope_row<Scalar, NX> const row = outputs.H.row(i) * power;
    Scalar const norm = row.norm();
    if(norm == Scalar{0})
        return std::nullopt;
    return std::pair{polytope_row<Scalar, NX>(row / norm), outputs.h(i) / norm};
}

template <typename Scalar, std::size_t NX>
void append_output_rows(polytopic_set<Scalar, NX>& set, const polytopic_set<Scalar, NX>& outputs, const Matrix<Scalar, NX, NX>& power)
{
    for(Eigen::Index i = 0; i < outputs.H.rows(); ++i)
        if(auto const unit = unit_output_row<Scalar, NX>(outputs, i, power))
            append_row<Scalar, NX>(set, unit->first, unit->second);
}

template <typename Scalar, std::size_t NX>
bool append_unimplied_rows(polytopic_set<Scalar, NX>& set, const polytopic_set<Scalar, NX>& outputs, const Matrix<Scalar, NX, NX>& power,
                           const std::vector<polytope_vertex<Scalar, NX>>& vertices, Scalar factor)
{
    bool appended = false;
    for(Eigen::Index i = 0; i < outputs.H.rows(); ++i)
    {
        auto const unit = unit_output_row<Scalar, NX>(outputs, i, power);
        if(!unit || implied<Scalar, NX>(unit->first, unit->second, vertices, factor))
            continue;
        append_row<Scalar, NX>(set, unit->first, unit->second);
        appended = true;
    }
    return appended;
}

/// Stein equation A_K^T X A_K - X = -I in Kronecker form, (A_K^T kron A_K^T - I) vec(X) = -vec(I). Its solution is
/// positive definite exactly when A_K is strictly stable, so it is the stability test as well as the Lyapunov matrix.
template <typename Scalar, std::size_t NX>
std::optional<Matrix<Scalar, NX, NX>> stein_lyapunov_matrix(const Matrix<Scalar, NX, NX>& A_K)
{
    constexpr int nx = static_cast<int>(NX);
    using kronecker = Eigen::Matrix<Scalar, nx * nx, nx * nx>;
    kronecker M;
    for(int i = 0; i < nx; ++i)
        for(int j = 0; j < nx; ++j)
            M.template block<nx, nx>(i * nx, j * nx) = A_K(j, i) * A_K.transpose();
    Eigen::FullPivLU<kronecker> const lu(M - kronecker::Identity());
    if(!lu.isInvertible())
        return std::nullopt;
    Matrix<Scalar, NX, NX> const identity = Matrix<Scalar, NX, NX>::Identity();
    Eigen::Matrix<Scalar, nx * nx, 1> const vec_X = lu.solve(-identity.reshaped());
    Matrix<Scalar, NX, NX> const X = vec_X.reshaped(nx, nx);
    if(!X.allFinite() || X.llt().info() != Eigen::Success)
        return std::nullopt;
    return X;
}

/// V(x) = x^T X x contracts by 1 - 1 / lambda_max(X) per step; with r2 the largest V on the vertices of O_{n-1}
/// and r1 the largest level of V inside Y, the last index O_inf needs is at most log(r1 / r2) / log(1 - 1 / lambda_max(X))
/// (@cite ossareh2023, Theorem 5, with the Stein matrix, so the bound exists for every strictly stable A_K). Rounded
/// up rather than down, so rounding in r1, r2 and lambda_max cannot cut it short.
template <typename Scalar, std::size_t NX>
Scalar admissible_step_bound(const Matrix<Scalar, NX, NX>& X, const polytopic_set<Scalar, NX>& outputs, const std::vector<polytope_vertex<Scalar, NX>>& vertices)
{
    auto const llt = X.llt();
    Scalar r1 = std::numeric_limits<Scalar>::infinity();
    for(Eigen::Index i = 0; i < outputs.H.rows(); ++i)
        r1 = std::min(r1, outputs.h(i) * outputs.h(i) / outputs.H.row(i).transpose().dot(llt.solve(outputs.H.row(i).transpose())));
    Scalar r2{0};
    for(auto const& v : vertices)
        r2 = std::max(r2, v.point.dot(X * v.point));
    if(r2 <= r1)
        return Scalar{0};
    Scalar const lambda_max = Eigen::SelfAdjointEigenSolver<Matrix<Scalar, NX, NX>>(X, Eigen::EigenvaluesOnly).eigenvalues().maxCoeff();
    return std::ceil(std::log(r1 / r2) / std::log1p(-Scalar{1} / lambda_max));
}

/// The first t whose rows C A_K^t are all implied ends the construction, since O_t = O_{t+1} implies O_t = O_inf.
template <typename Scalar, std::size_t NX>
ctrlpp::expected<polytopic_set<Scalar, NX>, terminal_set_error> grow_to_fixed_point(polytopic_set<Scalar, NX> set, const polytopic_set<Scalar, NX>& outputs,
                                                                                    const Matrix<Scalar, NX, NX>& A_K, Matrix<Scalar, NX, NX> power,
                                                                                    const Matrix<Scalar, NX, NX>& X, Scalar factor)
{
    auto vertices = enumerate_vertices<Scalar, NX>(set, factor);
    Scalar const cap = std::max(static_cast<Scalar>(NX), admissible_step_bound<Scalar, NX>(X, outputs, vertices));
    for(std::int64_t t = static_cast<std::int64_t>(NX);; ++t)
    {
        set = minimal_rows<Scalar, NX>(set, vertices, factor);
        if(!append_unimplied_rows<Scalar, NX>(set, outputs, power, vertices, factor))
            return set;
        if(static_cast<Scalar>(t) > cap)
            return ctrlpp::unexpected(terminal_set_error::not_converged);
        power = power * A_K;
        vertices = enumerate_vertices<Scalar, NX>(set, factor);
    }
}

}

/// @brief The maximal output admissible set O_inf = {x : C A_K^k x in Y for all k >= 0} of x+ = A_K x, with
/// C = [I; K] and Y = {x_min <= x <= x_max, u_min <= K x <= u_max} (@cite gilbert1991).
///
/// O_inf is positively invariant under A_K and the largest such set inside Y, so with A_K = A + B K it is the
/// terminal set of the LQR law u = K x. Infinite bound sides are absent. The rows C A_K^k up to k = NX - 1 bound
/// the set; each further k appends the rows its vertices do not imply (vertex enumeration, no linear program), and
/// redundant rows are pruned every step, so the result is a minimal H-representation with unit rows.
/// `rounding_factor` (positive, order one) scales the prefactor of every rounding bound. Rejections, in order:
///  * a finite bound side that does not strictly straddle zero, or a NaN bound -> origin_not_interior
///  * A_K not strictly stable, or A_K or K not finite                          -> closed_loop_not_stable
///  * the rows C A_K^k, k < NX, leave a direction in which the set recedes     -> unbounded_terminal_set
///  * no fixed point within the Lyapunov step bound                            -> not_converged
template <typename Scalar, std::size_t NX, std::size_t NU>
ctrlpp::expected<polytopic_set<Scalar, NX>, terminal_set_error> compute_output_admissible_set(
    const Matrix<Scalar, NX, NX>& A_K, const Matrix<Scalar, NU, NX>& K, const Vector<Scalar, NU>& u_min, const Vector<Scalar, NU>& u_max,
    const Vector<Scalar, NX>& x_min, const Vector<Scalar, NX>& x_max, Scalar rounding_factor = Scalar{1})
{
    static_assert(NX <= 4, "Output admissible set computation restricted to NX <= 4");
    if(!detail::straddles_zero<Scalar, NX>(x_min, x_max) || !detail::straddles_zero<Scalar, NU>(u_min, u_max))
        return ctrlpp::unexpected(terminal_set_error::origin_not_interior);
    auto const X = detail::stein_lyapunov_matrix<Scalar, NX>(A_K);
    if(!X || !K.allFinite())
        return ctrlpp::unexpected(terminal_set_error::closed_loop_not_stable);
    auto const outputs = detail::finite_outputs<Scalar, NX, NU>(K, u_min, u_max, x_min, x_max);
    polytopic_set<Scalar, NX> set{.H = Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>(0, static_cast<Eigen::Index>(NX)), .h = Eigen::VectorX<Scalar>(0)};
    Matrix<Scalar, NX, NX> power = Matrix<Scalar, NX, NX>::Identity();
    for(std::size_t k = 0; k < NX; ++k, power = power * A_K)
        detail::append_output_rows<Scalar, NX>(set, outputs, power);
    if(set.H.rows() < static_cast<Eigen::Index>(NX) || detail::has_recession_direction<Scalar, NX>(set, rounding_factor))
        return ctrlpp::unexpected(terminal_set_error::unbounded_terminal_set);
    return detail::grow_to_fixed_point<Scalar, NX>(std::move(set), outputs, A_K, power, *X, rounding_factor);
}

}

#endif
