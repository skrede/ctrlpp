#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_POLYTOPE_VERTICES_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_POLYTOPE_VERTICES_H

#include "ctrlpp/mpc/terminal_set.h"

#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <vector>
#include <cstddef>
#include <numeric>
#include <optional>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NU>
std::optional<Vector<Scalar, NU>> try_vertex_at_intersection(const polytopic_set<Scalar, NU>& constraints, const std::vector<Eigen::Index>& indices, Scalar tol)
{
    constexpr int nu = static_cast<int>(NU);
    Eigen::Index m = constraints.H.rows();

    Eigen::Matrix<Scalar, nu, nu> H_sub;
    Eigen::Matrix<Scalar, nu, 1> h_sub;
    for(Eigen::Index i = 0; i < nu; ++i)
    {
        H_sub.row(i) = constraints.H.row(indices[static_cast<std::size_t>(i)]);
        h_sub(i) = constraints.h(indices[static_cast<std::size_t>(i)]);
    }

    auto qr = H_sub.colPivHouseholderQr();
    if(!qr.isInvertible())
        return std::nullopt;

    Vector<Scalar, NU> v = qr.solve(h_sub);

    for(Eigen::Index j = 0; j < m; ++j)
        if(constraints.H.row(j).dot(v) > constraints.h(j) + tol)
            return std::nullopt;

    return v;
}

template <typename Scalar, std::size_t NU>
std::vector<Vector<Scalar, NU>> enumerate_polytope_vertices(const polytopic_set<Scalar, NU>& constraints, Scalar tol)
{
    constexpr auto nu = static_cast<Eigen::Index>(NU);
    std::vector<Vector<Scalar, NU>> vertices;
    Eigen::Index m = constraints.H.rows();
    if(m < nu)
        return vertices;
    std::vector<Eigen::Index> indices(NU);
    std::iota(indices.begin(), indices.end(), Eigen::Index{0});
    auto next_combination = [&]() -> bool
    {
        for(Eigen::Index i = nu - 1; i >= 0; --i)
            if(auto idx = static_cast<std::size_t>(i); indices[idx] < m - nu + i)
            {
                ++indices[idx];
                for(Eigen::Index j = i + 1; j < nu; ++j)
                    indices[static_cast<std::size_t>(j)] = indices[static_cast<std::size_t>(j - 1)] + 1;
                return true;
            }
        return false;
    };
    do
        if(auto v = try_vertex_at_intersection<Scalar, NU>(constraints, indices, tol))
            vertices.push_back(*v);
    while(next_combination());
    return vertices;
}

}

#endif
