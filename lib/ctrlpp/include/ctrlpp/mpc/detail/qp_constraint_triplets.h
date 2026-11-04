#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_QP_CONSTRAINT_TRIPLETS_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_QP_CONSTRAINT_TRIPLETS_H

#include "ctrlpp/mpc/terminal_set.h"

#include "ctrlpp/types.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>
#include <cstddef>
#include <cstdint>
#include <variant>
#include <type_traits>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX>
inline void add_initial_state_triplets(std::vector<Eigen::Triplet<Scalar>>& trips, Eigen::Index& row)
{
    constexpr int nx = static_cast<int>(NX);
    for(Eigen::Index i = 0; i < nx; ++i)
        trips.emplace_back(row + i, i, Scalar{1});
    row += nx;
}

template <typename Scalar, std::size_t NX, std::size_t NU>
inline void add_dynamics_triplets(
        std::vector<Eigen::Triplet<Scalar>>& trips, Eigen::Index& row, std::int32_t N, Eigen::Index n_x_total, const Matrix<Scalar, NX, NX>& A_sys, const Matrix<Scalar, NX, NU>& B_sys)
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);

    add_initial_state_triplets<Scalar, NX>(trips, row);
    for(std::int32_t k = 1; k <= N; ++k)
    {
        Eigen::Index x_k = k * nx;
        Eigen::Index x_km1 = (k - 1) * nx;
        Eigen::Index u_km1 = n_x_total + (k - 1) * nu;

        for(Eigen::Index i = 0; i < nx; ++i)
            for(Eigen::Index j = 0; j < nx; ++j)
                if(A_sys(i, j) != Scalar{0})
                    trips.emplace_back(row + i, x_km1 + j, -A_sys(i, j));
        for(Eigen::Index i = 0; i < nx; ++i)
            trips.emplace_back(row + i, x_k + i, Scalar{1});
        for(Eigen::Index i = 0; i < nx; ++i)
            for(Eigen::Index j = 0; j < nu; ++j)
                if(B_sys(i, j) != Scalar{0})
                    trips.emplace_back(row + i, u_km1 + j, -B_sys(i, j));
        row += nx;
    }
}

template <typename Scalar, std::size_t NX, std::size_t NU>
inline void add_state_bound_triplets(std::vector<Eigen::Triplet<Scalar>>& trips, Eigen::Index& row, std::int32_t N, Eigen::Index n_x_total, Eigen::Index n_u_total, bool has_soft)
{
    constexpr int nx = static_cast<int>(NX);
    for(std::int32_t k = 1; k <= N; ++k)
    {
        Eigen::Index x_k = k * nx;
        for(Eigen::Index i = 0; i < nx; ++i)
            trips.emplace_back(row + i, x_k + i, Scalar{1});

        if(has_soft)
        {
            Eigen::Index slack = n_x_total + n_u_total + (k - 1) * nx;
            for(Eigen::Index i = 0; i < nx; ++i)
                trips.emplace_back(row + i, slack + i, Scalar{-1});
        }
        row += nx;
    }
}

template <typename Scalar, std::size_t NX, std::size_t NU>
inline void add_input_bound_triplets(std::vector<Eigen::Triplet<Scalar>>& trips, Eigen::Index& row, std::int32_t N, Eigen::Index n_x_total)
{
    constexpr int nu = static_cast<int>(NU);
    for(std::int32_t k = 0; k < N; ++k)
    {
        Eigen::Index u_k = n_x_total + k * nu;
        for(Eigen::Index i = 0; i < nu; ++i)
            trips.emplace_back(row + i, u_k + i, Scalar{1});
        row += nu;
    }
}

template <typename Scalar, std::size_t NX, std::size_t NU>
inline void add_rate_bound_triplets(std::vector<Eigen::Triplet<Scalar>>& trips, Eigen::Index& row, std::int32_t N, Eigen::Index n_x_total)
{
    constexpr int nu = static_cast<int>(NU);
    for(std::int32_t k = 0; k < N; ++k)
    {
        Eigen::Index u_k = n_x_total + k * nu;
        for(Eigen::Index i = 0; i < nu; ++i)
            trips.emplace_back(row + i, u_k + i, Scalar{1});
        if(k > 0)
        {
            Eigen::Index u_km1 = n_x_total + (k - 1) * nu;
            for(Eigen::Index i = 0; i < nu; ++i)
                trips.emplace_back(row + i, u_km1 + i, Scalar{-1});
        }
        row += nu;
    }
}

template <typename Scalar, std::size_t NX>
inline void add_ellipsoidal_triplets(std::vector<Eigen::Triplet<Scalar>>& trips, Eigen::Index& row, Eigen::Index x_N_offset, const ellipsoidal_set<Scalar, NX>& s)
{
    constexpr int nx = static_cast<int>(NX);
    Eigen::SelfAdjointEigenSolver<Matrix<Scalar, NX, NX>> eig(s.P);
    auto V = eig.eigenvectors();

    // One two-sided row v_i^T x in [-bound, +bound] per eigendirection encodes the box; a row on -v_i^T x would
    // repeat it, since the two-sided box on -v_i^T x is the same set.
    for(Eigen::Index i = 0; i < nx; ++i)
        for(Eigen::Index j = 0; j < nx; ++j)
            if(V(j, i) != Scalar{0})
                trips.emplace_back(row + i, x_N_offset + j, V(j, i));
    row += nx;
}

template <typename Scalar, std::size_t NX>
inline void add_polytopic_triplets(std::vector<Eigen::Triplet<Scalar>>& trips, Eigen::Index& row, Eigen::Index x_N_offset, const polytopic_set<Scalar, NX>& s)
{
    constexpr int nx = static_cast<int>(NX);
    Eigen::Index n_faces = s.H.rows();
    for(Eigen::Index i = 0; i < n_faces; ++i)
        for(Eigen::Index j = 0; j < nx; ++j)
            if(s.H(i, j) != Scalar{0})
                trips.emplace_back(row + i, x_N_offset + j, s.H(i, j));
    row += n_faces;
}

template <typename Scalar, std::size_t NX>
inline void add_terminal_set_triplets(std::vector<Eigen::Triplet<Scalar>>& trips, Eigen::Index& row, Eigen::Index x_N_offset, const terminal_set<Scalar, NX>& tset)
{
    std::visit(
            [&](const auto& s)
            {
                using T = std::decay_t<decltype(s)>;
                if constexpr(std::is_same_v<T, ellipsoidal_set<Scalar, NX>>)
                    add_ellipsoidal_triplets<Scalar, NX>(trips, row, x_N_offset, s);
                else
                    add_polytopic_triplets<Scalar, NX>(trips, row, x_N_offset, s);
            },
            tset);
}

}

#endif
