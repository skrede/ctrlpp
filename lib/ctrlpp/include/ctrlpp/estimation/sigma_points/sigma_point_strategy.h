#ifndef HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_SIGMA_POINT_STRATEGY_H
#define HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_SIGMA_POINT_STRATEGY_H

#include "ctrlpp/types.h"

#include <array>
#include <cstddef>
#include <concepts>

namespace ctrlpp
{

template <typename Scalar, std::size_t NX, std::size_t NumPoints>
struct sigma_result
{
    std::array<Vector<Scalar, NX>, NumPoints> points;
    std::array<Scalar, NumPoints> Wm;
    std::array<Scalar, NumPoints> Wc;
};

template <typename S, typename Scalar, std::size_t NX>
concept sigma_point_strategy = requires {
    { S::num_points } -> std::convertible_to<std::size_t>;
    typename S::options_t;
} && requires(const S& s, const Vector<Scalar, NX>& x, const Matrix<Scalar, NX, NX>& P) {
    { s.generate(x, P) } -> std::convertible_to<sigma_result<Scalar, NX, S::num_points>>;
};

} // namespace ctrlpp

#endif
