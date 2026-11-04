#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_INVARIANT_SET_ROWS_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_INVARIANT_SET_ROWS_H

#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <vector>
#include <cstddef>
#include <utility>

namespace ctrlpp::detail
{

// The pre-image is stacked over every input vertex, { x : forall u in U, A x + B u in target }, so the iteration
// yields a robust control-invariant set, not the maximal (exists-u) one.
template <typename Scalar, std::size_t NX, std::size_t NU>
std::pair<Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>, Eigen::VectorX<Scalar>> compute_robust_pre_image(
        const Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>& H_curr, const Eigen::VectorX<Scalar>& h_curr, const Matrix<Scalar, NX, NX>& A_sys,
        const Matrix<Scalar, NX, NU>& B_sys, const std::vector<Vector<Scalar, NU>>& u_vertices)
{
    constexpr int nx = static_cast<int>(NX);
    auto H_A = (H_curr * A_sys).eval();
    Eigen::Index n_faces = H_curr.rows();
    auto n_verts = static_cast<Eigen::Index>(u_vertices.size());

    Eigen::Matrix<Scalar, Eigen::Dynamic, nx> H_pre(n_faces * n_verts, nx);
    Eigen::VectorX<Scalar> h_pre(n_faces * n_verts);

    for(Eigen::Index v = 0; v < n_verts; ++v)
    {
        auto Bu = (B_sys * u_vertices[static_cast<std::size_t>(v)]).eval();
        auto offset = (H_curr * Bu).eval();

        H_pre.middleRows(v * n_faces, n_faces) = H_A;
        h_pre.segment(v * n_faces, n_faces) = h_curr - offset;
    }

    return {H_pre, h_pre};
}

// Only a parallel halfplane no tighter than a kept one counts as redundant; no linear program is solved.
template <typename Scalar, std::size_t NX>
bool is_halfplane_redundant(const Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>& H, const Eigen::VectorX<Scalar>& h, Eigen::Index i, const std::vector<Eigen::Index>& keep,
        Scalar convergence_tol)
{
    Scalar row_norm = H.row(i).norm();
    Scalar h_normalized = h(i) / row_norm;
    Vector<Scalar, NX> dir_i = H.row(i).transpose() / row_norm;

    for(Eigen::Index j : keep)
    {
        Scalar row_norm_j = H.row(j).norm();
        Vector<Scalar, NX> dir_j = H.row(j).transpose() / row_norm_j;
        if(dir_i.dot(dir_j) > Scalar{1} - Scalar{1e-8})
        {
            Scalar h_norm_j = h(j) / row_norm_j;
            if(h_normalized >= h_norm_j - convergence_tol)
                return true;
        }
    }
    return false;
}

template <typename Scalar, std::size_t NX>
std::pair<Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>, Eigen::VectorX<Scalar>> extract_rows(
        const Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>& H, const Eigen::VectorX<Scalar>& h, const std::vector<Eigen::Index>& keep)
{
    constexpr int nx = static_cast<int>(NX);
    Eigen::Matrix<Scalar, Eigen::Dynamic, nx> H_out(static_cast<Eigen::Index>(keep.size()), nx);
    Eigen::VectorX<Scalar> h_out(static_cast<Eigen::Index>(keep.size()));
    for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(keep.size()); ++i)
    {
        H_out.row(i) = H.row(keep[static_cast<std::size_t>(i)]);
        h_out(i) = h(keep[static_cast<std::size_t>(i)]);
    }
    return {H_out, h_out};
}

// A truncated set is only an outer approximation of the filtered one.
template <typename Scalar, std::size_t NX>
struct filtered_halfplanes
{
    Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)> H;
    Eigen::VectorX<Scalar> h;
    bool truncated;
};

template <typename Scalar, std::size_t NX>
filtered_halfplanes<Scalar, NX> filter_redundant_halfplanes(
        const Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>& H, const Eigen::VectorX<Scalar>& h, Scalar convergence_tol)
{
    Eigen::Index n_rows = H.rows();

    std::vector<Eigen::Index> keep;
    keep.reserve(static_cast<std::size_t>(n_rows));

    for(Eigen::Index i = 0; i < n_rows; ++i)
    {
        if(H.row(i).norm() < Scalar{1e-14})
            continue;
        if(!is_halfplane_redundant<Scalar, NX>(H, h, i, keep, convergence_tol))
            keep.push_back(i);
    }

    // Bounded-memory guard on the halfplane count; the truncation is reported, never silent.
    constexpr std::size_t max_halfplanes = 500;
    bool truncated = keep.size() > max_halfplanes;
    if(truncated)
        keep.resize(max_halfplanes);

    auto [H_out, h_out] = extract_rows<Scalar, NX>(H, h, keep);
    return {.H = std::move(H_out), .h = std::move(h_out), .truncated = truncated};
}

}

#endif
