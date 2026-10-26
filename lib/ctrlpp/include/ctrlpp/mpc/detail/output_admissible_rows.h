#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_OUTPUT_ADMISSIBLE_ROWS_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_OUTPUT_ADMISSIBLE_ROWS_H

/// @brief Output rows C A_K^t, the step bound and the vertex certificate of the maximal output admissible set.
///
/// @cite gilbert1991 -- Gilbert & Tan, "Linear systems with state and control constraints: the theory and application of maximal output admissible sets", IEEE TAC 36(9), 1991
/// @cite ossareh2023 -- Ossareh & Kolmanovsky, "On complexity bounds for the maximal admissible set of linear time-invariant systems", 2023 (a-priori step bound)

#include "ctrlpp/mpc/detail/polytope_recession.h"
#include "ctrlpp/mpc/detail/polytope_redundancy.h"

#include "ctrlpp/mpc/terminal_set.h"

#include "ctrlpp/detail/covariance_ops.h"

#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <optional>
#include <algorithm>

namespace ctrlpp::detail
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

/// A_K v with its forward error: the vertex's own error mapped through |A_K|, plus each component's NX-term inner-product
/// bound (@cite higham2002, Sec. 3.1), which is what a row of A_K that cancels heavily contributes.
template <typename Scalar, std::size_t NX>
std::vector<polytope_vertex<Scalar, NX>> mapped_vertices(const std::vector<polytope_vertex<Scalar, NX>>& vertices, const Matrix<Scalar, NX, NX>& A_K, Scalar factor)
{
    std::vector<polytope_vertex<Scalar, NX>> images;
    images.reserve(vertices.size());
    for(auto const& v : vertices)
    {
        Vector<Scalar, NX> const magnitude = A_K.cwiseAbs() * v.point.cwiseAbs();
        Vector<Scalar, NX> error = A_K.cwiseAbs() * v.error;
        for(Eigen::Index k = 0; k < static_cast<Eigen::Index>(NX); ++k)
            error(k) += rounding_bound<Scalar>(static_cast<std::int32_t>(NX), magnitude(k), factor);
        images.push_back(polytope_vertex<Scalar, NX>{.point = A_K * v.point, .error = error});
    }
    return images;
}

/// The set, checked on its own vertices, recedes in no direction, lies inside Y, and contains the image under A_K of
/// every vertex, so A_K maps it into itself. Each implication carries the comparison bound the construction decides
/// with, so this checks the result against Y and A_K rather than adding a tolerance model. A bounded full-dimensional
/// polytope has at least NX + 1 vertices.
template <typename Scalar, std::size_t NX>
bool certified_admissible(const polytopic_set<Scalar, NX>& set, const polytopic_set<Scalar, NX>& outputs, const Matrix<Scalar, NX, NX>& A_K, Scalar factor)
{
    auto const vertices = enumerate_vertices<Scalar, NX>(set, factor);
    if(vertices.size() < NX + 1 || has_recession_direction<Scalar, NX>(set, factor))
        return false;
    Matrix<Scalar, NX, NX> const identity = Matrix<Scalar, NX, NX>::Identity();
    for(Eigen::Index i = 0; i < outputs.H.rows(); ++i)
    {
        auto const unit = unit_output_row<Scalar, NX>(outputs, i, identity);
        if(unit && !implied<Scalar, NX>(unit->first, unit->second, vertices, factor))
            return false;
    }
    auto const images = mapped_vertices<Scalar, NX>(vertices, A_K, factor);
    for(Eigen::Index i = 0; i < set.H.rows(); ++i)
        if(!implied<Scalar, NX>(set.H.row(i), set.h(i), images, factor))
            return false;
    return true;
}
}

#endif
