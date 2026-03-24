#ifndef HPP_GUARD_CTRLPP_CONTROL_LQR_H
#define HPP_GUARD_CTRLPP_CONTROL_LQR_H

/// @brief Linear Quadratic Regulator: infinite/finite/time-varying/integral-action.
///
/// @cite anderson1990 -- Anderson & Moore, "Optimal Control: Linear Quadratic Methods", 1990

#include "ctrlpp/control/dare.h"
#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <span>
#include <vector>
#include <cstddef>
#include <utility>
#include <optional>

namespace ctrlpp
{

// LQI gain result: partitioned feedback gain for integral action.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct lqi_result
{
    Eigen::Matrix<Scalar, int(NU), int(NX)> Kx;
    Eigen::Matrix<Scalar, int(NU), int(NY)> Ki;
};

// Infinite-horizon LQR gain via DARE.
// Returns K = (R + B^T P B)^{-1} B^T P A where P solves the DARE.
template <typename Scalar, std::size_t NX, std::size_t NU>
std::optional<Eigen::Matrix<Scalar, int(NU), int(NX)>>
lqr_gain(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B, const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
{
    auto P_opt = dare<Scalar, NX, NU>(A, B, Q, R);
    if(!P_opt)
        return std::nullopt;

    auto& P = *P_opt;
    auto BtP = (B.transpose() * P).eval();
    auto S = (R + BtP * B).eval();
    auto K = S.colPivHouseholderQr().solve(BtP * A).eval();
    return K;
}

// Infinite-horizon LQR gain with cross-weight N.
// K = (R + B^T P B)^{-1} (B^T P A + N^T)
template <typename Scalar, std::size_t NX, std::size_t NU>
std::optional<Eigen::Matrix<Scalar, int(NU), int(NX)>> lqr_gain(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
                                                                const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
                                                                const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
                                                                const Eigen::Matrix<Scalar, int(NU), int(NU)>& R,
                                                                const Eigen::Matrix<Scalar, int(NX), int(NU)>& N)
{
    auto P_opt = dare<Scalar, NX, NU>(A, B, Q, R, N);
    if(!P_opt)
        return std::nullopt;

    auto& P = *P_opt;
    auto BtP = (B.transpose() * P).eval();
    auto S = (R + BtP * B).eval();
    Eigen::Matrix<Scalar, int(NU), int(NX)> rhs = (BtP * A + N.transpose()).eval();
    auto K = S.colPivHouseholderQr().solve(rhs).eval();
    return K;
}

// Finite-horizon LQR via backward Riccati recursion.
// Returns gain sequence {K_0, K_1, ..., K_{N-1}} indexed by time step.
template <typename Scalar, std::size_t NX, std::size_t NU>
std::vector<Eigen::Matrix<Scalar, int(NU), int(NX)>> lqr_finite(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
                                                                const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
                                                                const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
                                                                const Eigen::Matrix<Scalar, int(NU), int(NU)>& R,
                                                                const Eigen::Matrix<Scalar, int(NX), int(NX)>& Qf,
                                                                std::size_t horizon)

{
    using MatNxN = Eigen::Matrix<Scalar, int(NX), int(NX)>;
    using MatK = Eigen::Matrix<Scalar, int(NU), int(NX)>;

    std::vector<MatK> gains(horizon);
    MatNxN P = Qf;

    for(std::size_t i = horizon; i > 0; --i)
    {
        auto BtP = (B.transpose() * P).eval();
        auto S = (R + BtP * B).eval();
        MatK K = S.colPivHouseholderQr().solve(BtP * A).eval();
        gains[i - 1] = K;
        P = (Q + A.transpose() * P * A - A.transpose() * P * B * K).eval();
    }

    return gains;
}

// Time-varying LQR via backward Riccati recursion with per-step matrices.
// Returns gain sequence {K_0, K_1, ..., K_{N-1}}.
template <typename Scalar, std::size_t NX, std::size_t NU>
std::vector<Eigen::Matrix<Scalar, int(NU), int(NX)>> lqr_tv_gains(const std::vector<Eigen::Matrix<Scalar, int(NX), int(NX)>>& As,
                                                                  const std::vector<Eigen::Matrix<Scalar, int(NX), int(NU)>>& Bs,
                                                                  const std::vector<Eigen::Matrix<Scalar, int(NX), int(NX)>>& Qs,
                                                                  const std::vector<Eigen::Matrix<Scalar, int(NU), int(NU)>>& Rs,
                                                                  const Eigen::Matrix<Scalar, int(NX), int(NX)>& Qf,
                                                                  std::size_t horizon)
{
    using MatNxN = Eigen::Matrix<Scalar, int(NX), int(NX)>;
    using MatK = Eigen::Matrix<Scalar, int(NU), int(NX)>;

    std::vector<MatK> gains(horizon);
    MatNxN P = Qf;

    for(std::size_t i = horizon; i > 0; --i)
    {
        std::size_t k = i - 1;
        const auto& Ak = As[k];
        const auto& Bk = Bs[k];
        const auto& Qk = Qs[k];
        const auto& Rk = Rs[k];

        auto BtP = (Bk.transpose() * P).eval();
        auto S = (Rk + BtP * Bk).eval();
        MatK K = S.colPivHouseholderQr().solve(BtP * Ak).eval();
        gains[k] = K;
        P = (Qk + Ak.transpose() * P * Ak - Ak.transpose() * P * Bk * K).eval();
    }

    return gains;
}

// LQI gain: augments state with integral of tracking error.
// Augmented system: A_aug = [[A, 0], [-C, I]], B_aug = [[B], [0]]
// Returns lqi_result with partitioned Kx (NU x NX) and Ki (NU x NY).
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
std::optional<lqi_result<Scalar, NX, NU, NY>> lqi_gain(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
                                                       const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
                                                       const Eigen::Matrix<Scalar, int(NY), int(NX)>& C,
                                                       const Eigen::Matrix<Scalar, int(NX + NY), int(NX + NY)>& Q_aug,
                                                       const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
{
    constexpr std::size_t NX_AUG = NX + NY;
    constexpr int nx = static_cast<int>(NX);
    constexpr int ny = static_cast<int>(NY);
    constexpr int nu = static_cast<int>(NU);
    constexpr int nx_aug = static_cast<int>(NX_AUG);

    using MatAug = Eigen::Matrix<Scalar, nx_aug, nx_aug>;
    using MatBaug = Eigen::Matrix<Scalar, nx_aug, nu>;

    // Build augmented A: [[A, 0], [-C, I]]
    MatAug A_aug = MatAug::Zero();
    A_aug.template block<nx, nx>(0, 0) = A;
    A_aug.template block<ny, nx>(nx, 0) = -C;
    A_aug.template block<ny, ny>(nx, nx) = Eigen::Matrix<Scalar, ny, ny>::Identity();

    // Build augmented B: [[B], [0]]
    MatBaug B_aug = MatBaug::Zero();
    B_aug.template block<nx, nu>(0, 0) = B;

    // Solve augmented LQR
    auto K_aug_opt = lqr_gain<Scalar, NX_AUG, NU>(A_aug, B_aug, Q_aug, R);
    if(!K_aug_opt)
        return std::nullopt;

    auto& K_aug = *K_aug_opt;

    lqi_result<Scalar, NX, NU, NY> result;
    result.Kx = K_aug.template block<nu, nx>(0, 0);
    result.Ki = K_aug.template block<nu, ny>(0, nx);
    return result;
}

// Evaluate quadratic trajectory cost: sum of x^T Q x + u^T R u.
// If xs has one more element than us, the terminal state cost x_N^T Q x_N is included.
template <typename Scalar, std::size_t NX, std::size_t NU>
auto lqr_cost(std::span<const Eigen::Matrix<Scalar, int(NX), 1>> xs,
              std::span<const Eigen::Matrix<Scalar, int(NU), 1>> us,
              const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
              const Eigen::Matrix<Scalar, int(NU), int(NU)>& R) -> Scalar
{
    Scalar cost{0};

    for(std::size_t k = 0; k < us.size(); ++k)
    {
        cost += (xs[k].transpose() * Q * xs[k])(0, 0);
        cost += (us[k].transpose() * R * us[k])(0, 0);
    }

    // Terminal state cost if xs has one more element
    if(xs.size() == us.size() + 1)
        cost += (xs.back().transpose() * Q * xs.back())(0, 0);

    return cost;
}

// Thin LQR controller class storing a precomputed gain matrix.
template <typename Scalar, std::size_t NX, std::size_t NU>
class lqr
{
public:
    using gain_type = Eigen::Matrix<Scalar, int(NU), int(NX)>;
    using state_type = Eigen::Matrix<Scalar, int(NX), 1>;
    using input_type = Eigen::Matrix<Scalar, int(NU), 1>;

    explicit lqr(gain_type K) : K_(std::move(K)) {}

    auto compute(const state_type& x) const -> input_type { return (-K_ * x).eval(); }

    auto gain() const -> const gain_type& { return K_; }

private:
    gain_type K_;
};

// Time-varying LQR controller storing a gain sequence.
template <typename Scalar, std::size_t NX, std::size_t NU>
class lqr_time_varying
{
public:
    using gain_type = Eigen::Matrix<Scalar, int(NU), int(NX)>;
    using state_type = Eigen::Matrix<Scalar, int(NX), 1>;
    using input_type = Eigen::Matrix<Scalar, int(NU), 1>;

    explicit lqr_time_varying(std::vector<gain_type> gains) : gains_(std::move(gains)) {}

    auto compute(const state_type& x, std::size_t k) const -> input_type { return (-gains_[k] * x).eval(); }

    auto gain(std::size_t k) const -> const gain_type& { return gains_[k]; }

    auto horizon() const -> std::size_t { return gains_.size(); }

private:
    std::vector<gain_type> gains_;
};

} // namespace ctrlpp

#endif
