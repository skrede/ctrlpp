#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_POLYTOPE_REDUNDANCY_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_POLYTOPE_REDUNDANCY_H

// In a bounded polytope {x : G x <= g} with unit rows, each comparison of a row with a vertex carries its own rounding
// bound: the residual G_i v - g_i as an (n + 1)-term inner product (@cite higham2002, Sec. 3.1), plus the vertex's
// forward error bounded from the residual of its n defining rows (@cite higham2002, Sec. 7.2). A looser feasibility
// bound only admits near-vertices, which keeps more rows and stays sound; the implication bound is what can admit a set
// larger than the exact one, by at most that bound.

#include "ctrlpp/mpc/terminal_set.h"

#include "ctrlpp/detail/covariance_ops.h"

#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <optional>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX>
using polytope_row = Eigen::Matrix<Scalar, 1, static_cast<int>(NX)>;

template <typename Scalar, std::size_t NX>
struct polytope_vertex
{
    Vector<Scalar, NX> point;
    Vector<Scalar, NX> error;
};

template <typename Scalar, std::size_t NX>
void append_row(polytopic_set<Scalar, NX>& set, const polytope_row<Scalar, NX>& row, Scalar bound)
{
    auto const m = set.H.rows();
    set.H.conservativeResize(m + 1, Eigen::NoChange);
    set.h.conservativeResize(m + 1);
    set.H.row(m) = row;
    set.h(m) = bound;
}

template <typename Scalar, std::size_t NX>
Scalar residual_bound(const polytope_row<Scalar, NX>& row, Scalar bound, const Vector<Scalar, NX>& point, Scalar factor)
{
    return rounding_bound<Scalar>(static_cast<std::int32_t>(NX + 1), row.cwiseAbs().dot(point.cwiseAbs()) + std::abs(bound), factor);
}

template <typename Scalar, std::size_t NX>
Scalar comparison_bound(const polytope_row<Scalar, NX>& row, Scalar bound, const polytope_vertex<Scalar, NX>& v, Scalar factor)
{
    return residual_bound<Scalar, NX>(row, bound, v.point, factor) + row.cwiseAbs().dot(v.error);
}

template <typename Scalar, std::size_t NX>
std::optional<polytope_vertex<Scalar, NX>> vertex_at(const polytopic_set<Scalar, NX>& set, const std::array<Eigen::Index, NX>& subset, Scalar factor)
{
    constexpr auto nx = static_cast<Eigen::Index>(NX);
    Matrix<Scalar, NX, NX> G;
    Vector<Scalar, NX> g;
    for(Eigen::Index k = 0; k < nx; ++k)
    {
        G.row(k) = set.H.row(subset[static_cast<std::size_t>(k)]);
        g(k) = set.h(subset[static_cast<std::size_t>(k)]);
    }
    Eigen::ColPivHouseholderQR<Matrix<Scalar, NX, NX>> const qr(G);
    if(!qr.isInvertible())
        return std::nullopt;
    Vector<Scalar, NX> const point = qr.solve(g);
    Vector<Scalar, NX> slack;
    for(Eigen::Index k = 0; k < nx; ++k)
        slack(k) = std::abs(G.row(k).dot(point) - g(k)) + residual_bound<Scalar, NX>(G.row(k), g(k), point, factor);
    return polytope_vertex<Scalar, NX>{.point = point, .error = qr.inverse().cwiseAbs() * slack};
}

template <typename Scalar, std::size_t NX>
bool admits(const polytopic_set<Scalar, NX>& set, const std::array<Eigen::Index, NX>& subset, const polytope_vertex<Scalar, NX>& v, Scalar factor)
{
    for(Eigen::Index j = 0; j < set.H.rows(); ++j)
        if(std::find(subset.begin(), subset.end(), j) == subset.end() && set.H.row(j).dot(v.point) - set.h(j) > comparison_bound<Scalar, NX>(set.H.row(j), set.h(j), v, factor))
            return false;
    return true;
}

template <std::size_t NX>
bool next_subset(std::array<Eigen::Index, NX>& subset, Eigen::Index m)
{
    for(std::size_t i = NX; i-- > 0;)
        if(subset[i] < m - static_cast<Eigen::Index>(NX - i))
        {
            ++subset[i];
            for(std::size_t j = i + 1; j < NX; ++j)
                subset[j] = subset[j - 1] + 1;
            return true;
        }
    return false;
}

// Vertex enumeration over every NX-row subset, so the cost is C(m, NX) small solves.
template <typename Scalar, std::size_t NX>
std::vector<polytope_vertex<Scalar, NX>> enumerate_vertices(const polytopic_set<Scalar, NX>& set, Scalar factor)
{
    std::vector<polytope_vertex<Scalar, NX>> vertices;
    if(set.H.rows() < static_cast<Eigen::Index>(NX))
        return vertices;
    std::array<Eigen::Index, NX> subset{};
    std::iota(subset.begin(), subset.end(), Eigen::Index{0});
    do
    {
        if(auto v = vertex_at<Scalar, NX>(set, subset, factor); v && admits<Scalar, NX>(set, subset, *v, factor))
            vertices.push_back(*v);
    } while(next_subset<NX>(subset, set.H.rows()));
    return vertices;
}

template <typename Scalar, std::size_t NX>
bool implied(const polytope_row<Scalar, NX>& row, Scalar bound, const std::vector<polytope_vertex<Scalar, NX>>& vertices, Scalar factor)
{
    return std::all_of(
            vertices.begin(), vertices.end(), [&](const polytope_vertex<Scalar, NX>& v) { return row.dot(v.point) - bound <= comparison_bound<Scalar, NX>(row, bound, v, factor); });
}

template <typename Scalar, std::size_t NX>
std::vector<std::size_t> tight_vertices(const polytope_row<Scalar, NX>& row, Scalar bound, const std::vector<polytope_vertex<Scalar, NX>>& vertices, Scalar factor)
{
    std::vector<std::size_t> tight;
    for(std::size_t k = 0; k < vertices.size(); ++k)
        if(row.dot(vertices[k].point) - bound >= -comparison_bound<Scalar, NX>(row, bound, vertices[k], factor))
            tight.push_back(k);
    return tight;
}

template <typename Scalar, std::size_t NX>
bool certifies_facet(
        const polytope_row<Scalar, NX>& row, const std::array<Eigen::Index, NX>& chosen, const std::vector<std::size_t>& tight, const std::vector<polytope_vertex<Scalar, NX>>& vertices)
{
    if constexpr(NX == 1)
        return true;
    else
    {
        constexpr int nx = static_cast<int>(NX);
        auto const& base = vertices[tight[static_cast<std::size_t>(chosen[0])]];
        Eigen::Matrix<Scalar, nx - 1, nx> projected;
        Scalar noise{0};
        for(std::size_t k = 1; k < NX; ++k)
        {
            auto const& v = vertices[tight[static_cast<std::size_t>(chosen[k])]];
            polytope_row<Scalar, NX> const offset = (v.point - base.point).transpose();
            projected.row(static_cast<Eigen::Index>(k - 1)) = offset - row.dot(offset) * row;
            noise = std::max(noise, v.error.norm() + base.error.norm());
        }
        return Eigen::JacobiSVD<Eigen::Matrix<Scalar, nx - 1, nx>>(projected).singularValues()(nx - 2) > std::sqrt(static_cast<Scalar>(nx - 1)) * noise;
    }
}

// A row is a facet exactly when NX of its tight vertices are affinely independent in its hyperplane: for each NX of them the offsets are
// projected onto it and sigma_{NX-1} is compared with the Frobenius bound of their forward errors (@cite golub2013, Sec. 8.6.1), so no badly located vertex can veto.
template <typename Scalar, std::size_t NX>
bool spans_facet(const polytope_row<Scalar, NX>& row, const std::vector<std::size_t>& tight, const std::vector<polytope_vertex<Scalar, NX>>& vertices)
{
    if(tight.size() < NX)
        return false;
    std::array<Eigen::Index, NX> chosen{};
    std::iota(chosen.begin(), chosen.end(), Eigen::Index{0});
    do
    {
        if(certifies_facet<Scalar, NX>(row, chosen, tight, vertices))
            return true;
    } while(next_subset<NX>(chosen, static_cast<Eigen::Index>(tight.size())));
    return false;
}

// Keeps the rows of `set` that are facets on its `vertices`; rows defining the same facet have the same tight vertices, so only the first is kept.
template <typename Scalar, std::size_t NX>
polytopic_set<Scalar, NX> minimal_rows(const polytopic_set<Scalar, NX>& set, const std::vector<polytope_vertex<Scalar, NX>>& vertices, Scalar factor)
{
    polytopic_set<Scalar, NX> kept{.H = Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>(0, static_cast<Eigen::Index>(NX)), .h = Eigen::VectorX<Scalar>(0)};
    std::vector<std::vector<std::size_t>> facets;
    for(Eigen::Index i = 0; i < set.H.rows(); ++i)
    {
        auto tight = tight_vertices<Scalar, NX>(set.H.row(i), set.h(i), vertices, factor);
        if(!spans_facet<Scalar, NX>(set.H.row(i), tight, vertices) || std::find(facets.begin(), facets.end(), tight) != facets.end())
            continue;
        facets.push_back(std::move(tight));
        append_row<Scalar, NX>(kept, set.H.row(i), set.h(i));
    }
    return kept;
}

}

#endif
