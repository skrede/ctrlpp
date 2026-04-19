#ifndef HPP_GUARD_CTRLPP_CONTROL_LQR_H
#define HPP_GUARD_CTRLPP_CONTROL_LQR_H

/// @brief Linear Quadratic Regulator: infinite/finite/time-varying/integral-action.
///
/// @cite anderson1990 -- Anderson & Moore, "Optimal Control: Linear Quadratic Methods", 1990

#include "ctrlpp/control/dare.h"
#include "ctrlpp/control/care.h"
#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <span>
#include <vector>
#include <cstddef>
#include <utility>
#include <optional>
#include <type_traits>

namespace ctrlpp
{

// LQI gain result: partitioned feedback gain for integral action.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct lqi_result
{
    static_assert(std::is_floating_point_v<Scalar>, "Scalar must be a floating-point type");
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NY > 0, "Output dimension NY must be positive");
    Eigen::Matrix<Scalar, int(NU), int(NX)> Kx;
    Eigen::Matrix<Scalar, int(NU), int(NY)> Ki;
};

// Infinite-horizon LQR gain via DARE.
// Returns K = (R + B^T P B)^{-1} B^T P A where P solves the DARE.
template <typename Scalar, std::size_t NX, std::size_t NU>
std::optional<Eigen::Matrix<Scalar, int(NU), int(NX)>>
lqr_gain(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B, const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
{
    auto P_result = dare<Scalar, NX, NU>(A, B, Q, R);
    if(!P_result)
        return std::nullopt;

    const auto& P = P_result->P;
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
    auto P_result = dare<Scalar, NX, NU>(A, B, Q, R, N);
    if(!P_result)
        return std::nullopt;

    const auto& P = P_result->P;
    auto BtP = (B.transpose() * P).eval();
    auto S = (R + BtP * B).eval();
    Eigen::Matrix<Scalar, int(NU), int(NX)> rhs = (BtP * A + N.transpose()).eval();
    auto K = S.colPivHouseholderQr().solve(rhs).eval();
    return K;
}

// Continuous-time infinite-horizon LQR gain via CARE.
// K = R^{-1} B^T P where P solves A^T P + P A - P B R^{-1} B^T P + Q = 0.
//
// Computes R^{-1} once via ldlt (R is SPD for valid LQR problems), builds the
// Hamiltonian using that pre-computed R^{-1}, and reuses it for the K formula --
// one matrix factorisation of R instead of two.
template <typename Scalar, std::size_t NX, std::size_t NU,
          detail::care_solve_method Method = detail::sign_function_care_method>
std::optional<Eigen::Matrix<Scalar, int(NU), int(NX)>>
lqr_gain_continuous(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
                    const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
                    const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
                    const Eigen::Matrix<Scalar, int(NU), int(NU)>& R,
                    Method                                         /*method_tag*/ = {})
{
    static_assert(std::is_floating_point_v<Scalar>, "Scalar must be a floating-point type");
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);
    constexpr int n2 = 2 * nx;
    using MatU    = Eigen::Matrix<Scalar, nu, nu>;
    using Mat2N   = Eigen::Matrix<Scalar, n2, n2>;

    if(!A.allFinite() || !B.allFinite() || !Q.allFinite() || !R.allFinite())
        return std::nullopt;

    const MatU R_inv = R.ldlt().solve(MatU::Identity());

    Mat2N H;
    H.template block<nx, nx>(0, 0) = A;
    H.template block<nx, nx>(0, nx) = -B * R_inv * B.transpose();
    H.template block<nx, nx>(nx, 0) = -Q;
    H.template block<nx, nx>(nx, nx) = -A.transpose();

    if(!H.allFinite())
        return std::nullopt;

    auto result = detail::care_solve_from_hamiltonian<Scalar, NX, Method>(H);
    if(!result)
        return std::nullopt;

    return (R_inv * B.transpose() * result->P).eval();
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

namespace detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
auto build_lqi_augmented_system(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
                                const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
                                const Eigen::Matrix<Scalar, int(NY), int(NX)>& C)
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
auto partition_lqi_gain(const Eigen::Matrix<Scalar, int(NX + NY), int(NX + NY)>& A_aug,
                        const Eigen::Matrix<Scalar, int(NX + NY), int(NU)>& B_aug,
                        const Eigen::Matrix<Scalar, int(NX + NY), int(NX + NY)>& Q_aug,
                        const Eigen::Matrix<Scalar, int(NU), int(NU)>& R) -> std::optional<lqi_result<Scalar, NX, NU, NY>>
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int ny = static_cast<int>(NY);
    constexpr int nu = static_cast<int>(NU);

    auto K_aug_opt = lqr_gain<Scalar, NX + NY, NU>(A_aug, B_aug, Q_aug, R);
    if(!K_aug_opt)
        return std::nullopt;

    auto& K_aug = *K_aug_opt;
    lqi_result<Scalar, NX, NU, NY> result;
    result.Kx = K_aug.template block<nu, nx>(0, 0);
    result.Ki = K_aug.template block<nu, ny>(0, nx);
    return result;
}

}

// LQI gain: augments state with integral of tracking error.
// Augmented system: A_aug = [[A, 0], [-C, I]], B_aug = [[B], [0]]
// Returns lqi_result with partitioned Kx (NU x NX) and Ki (NU x NY).
/// @cite anderson1990 -- Anderson & Moore, "Optimal Control: Linear Quadratic Methods", 1990, Ch. 9 (integral action)
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
std::optional<lqi_result<Scalar, NX, NU, NY>> lqi_gain(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
                                                       const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
                                                       const Eigen::Matrix<Scalar, int(NY), int(NX)>& C,
                                                       const Eigen::Matrix<Scalar, int(NX + NY), int(NX + NY)>& Q_aug,
                                                       const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
{
    auto [A_aug, B_aug] = detail::build_lqi_augmented_system<Scalar, NX, NU, NY>(A, B, C);
    return detail::partition_lqi_gain<Scalar, NX, NU, NY>(A_aug, B_aug, Q_aug, R);
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
    static_assert(std::is_floating_point_v<Scalar>, "Scalar must be a floating-point type");
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

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
    static_assert(std::is_floating_point_v<Scalar>, "Scalar must be a floating-point type");
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

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

}

#endif
