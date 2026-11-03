#ifndef HPP_GUARD_CTRLPP_CONTROL_DETAIL_LQR_DISCRETE_H
#define HPP_GUARD_CTRLPP_CONTROL_DETAIL_LQR_DISCRETE_H

#include "ctrlpp/control/lqr_types.h"

#include "ctrlpp/expected.h"

#include <Eigen/Dense>

#include <span>
#include <vector>
#include <cstddef>
#include <utility>
#include <algorithm>

namespace ctrlpp::detail
{

// One step of the backward Riccati recursion: returns K and replaces P with its predecessor.
template <typename Scalar, std::size_t NX, std::size_t NU>
Eigen::Matrix<Scalar, int(NU), int(NX)> lqr_riccati_step(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R, Eigen::Matrix<Scalar, int(NX), int(NX)>& P)
{
    auto BtP = (B.transpose() * P).eval();
    auto S = (R + BtP * B).eval();
    Eigen::Matrix<Scalar, int(NU), int(NX)> K = S.colPivHouseholderQr().solve(BtP * A).eval();
    P = (Q + A.transpose() * P * A - A.transpose() * P * B * K).eval();
    return K;
}

template <typename Scalar, std::size_t NX, std::size_t NU>
ctrlpp::expected<void, lqr_error> check_lqr_tv_sequences(const std::vector<Eigen::Matrix<Scalar, int(NX), int(NX)>>& As, const std::vector<Eigen::Matrix<Scalar, int(NX), int(NU)>>& Bs,
        const std::vector<Eigen::Matrix<Scalar, int(NX), int(NX)>>& Qs, const std::vector<Eigen::Matrix<Scalar, int(NU), int(NU)>>& Rs,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& Qf)
{
    if(Bs.size() != As.size() || Qs.size() != As.size() || Rs.size() != As.size())
        return ctrlpp::unexpected(lqr_error::sequence_length_mismatch);

    const auto finite = [](const auto& M) { return M.allFinite(); };
    if(!Qf.allFinite() || !std::ranges::all_of(As, finite) || !std::ranges::all_of(Bs, finite) || !std::ranges::all_of(Qs, finite) || !std::ranges::all_of(Rs, finite))
        return ctrlpp::unexpected(lqr_error::non_finite_input);
    return {};
}

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
std::pair<Eigen::Matrix<Scalar, int(NX + NY), int(NX + NY)>, Eigen::Matrix<Scalar, int(NX + NY), int(NU)>> build_lqi_augmented_system(
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B, const Eigen::Matrix<Scalar, int(NY), int(NX)>& C)
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int ny = static_cast<int>(NY);
    constexpr int nu = static_cast<int>(NU);
    constexpr int nx_aug = static_cast<int>(NX + NY);

    using MatAug = Eigen::Matrix<Scalar, nx_aug, nx_aug>;
    using MatBaug = Eigen::Matrix<Scalar, nx_aug, nu>;

    MatAug A_aug = MatAug::Zero();
    A_aug.template block<nx, nx>(0, 0) = A;
    A_aug.template block<ny, nx>(nx, 0) = -C;
    A_aug.template block<ny, ny>(nx, nx) = Eigen::Matrix<Scalar, ny, ny>::Identity();

    MatBaug B_aug = MatBaug::Zero();
    B_aug.template block<nx, nu>(0, 0) = B;

    return std::pair{A_aug, B_aug};
}

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
lqi_result<Scalar, NX, NU, NY> partition_lqi_gain(const Eigen::Matrix<Scalar, int(NU), int(NX + NY)>& K_aug)
{
    lqi_result<Scalar, NX, NU, NY> result;
    result.Kx = K_aug.template block<int(NU), int(NX)>(0, 0);
    result.Ki = K_aug.template block<int(NU), int(NY)>(0, int(NX));
    return result;
}

template <typename Scalar, std::size_t NX, std::size_t NU>
bool lqr_cost_operands_finite(std::span<const Eigen::Matrix<Scalar, int(NX), 1>> xs, std::span<const Eigen::Matrix<Scalar, int(NU), 1>> us,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
{
    const auto finite = [](const auto& v) { return v.allFinite(); };
    return Q.allFinite() && R.allFinite() && std::ranges::all_of(xs, finite) && std::ranges::all_of(us, finite);
}

// A terminal x_N' Q x_N is added when xs holds one more state than us.
template <typename Scalar, std::size_t NX, std::size_t NU>
Scalar lqr_trajectory_cost(std::span<const Eigen::Matrix<Scalar, int(NX), 1>> xs, std::span<const Eigen::Matrix<Scalar, int(NU), 1>> us,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
{
    Scalar cost{0};
    for(std::size_t k = 0; k < us.size(); ++k)
    {
        cost += (xs[k].transpose() * Q * xs[k])(0, 0);
        cost += (us[k].transpose() * R * us[k])(0, 0);
    }

    if(xs.size() == us.size() + 1)
        cost += (xs.back().transpose() * Q * xs.back())(0, 0);
    return cost;
}

}

#endif
