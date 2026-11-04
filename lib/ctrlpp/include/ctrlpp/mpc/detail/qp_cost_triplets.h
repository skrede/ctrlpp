#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_QP_COST_TRIPLETS_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_QP_COST_TRIPLETS_H

#include "ctrlpp/types.h"

#include <Eigen/Sparse>

#include <vector>
#include <cstddef>
#include <cstdint>
#include <optional>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU>
inline void add_diagonal_block(std::vector<Eigen::Triplet<Scalar>>& trips, Eigen::Index offset, const auto& block, Eigen::Index size)
{
    for(Eigen::Index i = 0; i < size; ++i)
        for(Eigen::Index j = i; j < size; ++j)
            if(block(i, j) != Scalar{0})
                trips.emplace_back(static_cast<int>(offset + i), static_cast<int>(offset + j), block(i, j));
}

template <typename Scalar, std::size_t NX, std::size_t NU>
inline void add_soft_constraint_triplets(
        std::vector<Eigen::Triplet<Scalar>>& trips, std::int32_t N, Eigen::Index slack_offset, Scalar soft_penalty, const std::optional<Vector<Scalar, NX>>& soft_state_penalty)
{
    constexpr int nx = static_cast<int>(NX);
    for(std::int32_t k = 0; k < N; ++k)
        for(Eigen::Index i = 0; i < nx; ++i)
        {
            Scalar penalty = soft_state_penalty.has_value() ? (*soft_state_penalty)(i) : soft_penalty;
            trips.emplace_back(static_cast<int>(slack_offset + k * nx + i), static_cast<int>(slack_offset + k * nx + i), penalty);
        }
}

}

#endif
