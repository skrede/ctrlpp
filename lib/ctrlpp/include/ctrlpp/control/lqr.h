#ifndef HPP_GUARD_CTRLPP_CONTROL_LQR_H
#define HPP_GUARD_CTRLPP_CONTROL_LQR_H

/// @brief Linear Quadratic Regulator: infinite/finite/time-varying/integral-action.
///
/// @cite anderson1990 -- Anderson & Moore, "Optimal Control: Linear Quadratic Methods", 1990

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/control/care.h"
#include "ctrlpp/control/dare.h"

#include <Eigen/Dense>

#include <span>
#include <limits>
#include <vector>
#include <cstddef>
#include <utility>

namespace ctrlpp
{

// LQI gain result: partitioned feedback gain for integral action.
template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct lqi_result
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NY > 0, "Output dimension NY must be positive");
    Eigen::Matrix<Scalar, int(NU), int(NX)> Kx;
    Eigen::Matrix<Scalar, int(NU), int(NY)> Ki;
};

/// Infinite-horizon LQR gain via DARE.
///
/// Returns K = (R + B^T P B)^{-1} B^T P A where P solves the DARE, or the
/// Riccati solver's own `dare_error` verbatim. The gain is a function of that
/// solve and has no failure mode of its own, so it forwards the enumerator
/// rather than restating the cause under a second name: an unstabilizable pair,
/// a singular state matrix and a non-converged factorization send the caller to
/// fix three different things.
template <typename Scalar, std::size_t NX, std::size_t NU>
auto lqr_gain(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B, const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
    -> ctrlpp::expected<Eigen::Matrix<Scalar, int(NU), int(NX)>, dare_error>
{
    auto P_result = dare<Scalar, NX, NU>(A, B, Q, R);
    if(!P_result)
        return ctrlpp::unexpected(P_result.error());

    const auto& P = P_result->P;
    auto BtP = (B.transpose() * P).eval();
    auto S = (R + BtP * B).eval();
    auto K = S.colPivHouseholderQr().solve(BtP * A).eval();
    return K;
}

/// Infinite-horizon LQR gain with cross-weight N.
///
/// K = (R + B^T P B)^{-1} (B^T P A + N^T). Forwards the Riccati solver's
/// `dare_error` for the same reason the cross-weight-free overload does.
template <typename Scalar, std::size_t NX, std::size_t NU>
auto lqr_gain(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
              const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
              const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
              const Eigen::Matrix<Scalar, int(NU), int(NU)>& R,
              const Eigen::Matrix<Scalar, int(NX), int(NU)>& N)
    -> ctrlpp::expected<Eigen::Matrix<Scalar, int(NU), int(NX)>, dare_error>
{
    auto P_result = dare<Scalar, NX, NU>(A, B, Q, R, N);
    if(!P_result)
        return ctrlpp::unexpected(P_result.error());

    const auto& P = P_result->P;
    auto BtP = (B.transpose() * P).eval();
    auto S = (R + BtP * B).eval();
    Eigen::Matrix<Scalar, int(NU), int(NX)> rhs = (BtP * A + N.transpose()).eval();
    auto K = S.colPivHouseholderQr().solve(rhs).eval();
    return K;
}

/// Continuous-time infinite-horizon LQR gain via CARE.
///
/// K = R^{-1} B^T P where P solves A^T P + P A - P B R^{-1} B^T P + Q = 0.
///
/// Computes R^{-1} once via ldlt (R is SPD for valid LQR problems), builds the
/// Hamiltonian using that pre-computed R^{-1}, and reuses it for the K formula --
/// one matrix factorisation of R instead of two.
///
/// Reports through `care_error`. It makes three rejections of its own before
/// calling anything: a non-finite argument and a Hamiltonian that overflowed
/// while being assembled are both `care_error::non_finite_input`, and a
/// rank-deficient R is `care_error::singular_r`. Everything else is the
/// solver's own enumerator forwarded.
template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU,
          detail::care_solve_method Method = detail::sign_function_care_method>
auto lqr_gain_continuous(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
                         const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
                         const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
                         const Eigen::Matrix<Scalar, int(NU), int(NU)>& R,
                         Method                                         /*method_tag*/ = {})
    -> ctrlpp::expected<Eigen::Matrix<Scalar, int(NU), int(NX)>, care_error>
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);
    constexpr int n2 = 2 * nx;
    using MatU    = Eigen::Matrix<Scalar, nu, nu>;
    using Mat2N   = Eigen::Matrix<Scalar, n2, n2>;

    if(!A.allFinite() || !B.allFinite() || !Q.allFinite() || !R.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    // This function forms R^{-1} itself rather than going through the Hamiltonian
    // build, so it needs its own singularity test -- and it needs one MORE than the
    // build does, because the factorization it uses fails quietly. Eigen's LDLT solve
    // zeroes the rank-deficient directions instead of producing infinities, so a
    // singular R yields a FINITE R^{-1} of zeros, an entirely finite Hamiltonian
    // describing a plant with no control authority at all, and the sign-function
    // iteration then stagnates on it. Nothing anywhere on that path is non-finite, so
    // no downstream check can catch it: measured, a zero R reported
    // sign_function_stagnated, which sends the caller to look at convergence when the
    // obstacle is the weighting they passed.
    //
    // The condition is read off the factorization already being formed. LDLT's pivots
    // are its diagonal, so the matrix is rank-deficient exactly when the smallest pivot
    // magnitude falls under the largest one scaled by the input dimension times unit
    // roundoff -- the same reciprocal-pivot convention, and the same size factor, the
    // two Riccati builds apply to their own inverted operands. Sign is deliberately not
    // tested: an indefinite but nonsingular R is a different condition from a singular
    // one and is left to whatever the solve path already does with it.
    const auto R_ldlt = R.ldlt();
    const auto R_pivots = R_ldlt.vectorD().cwiseAbs().eval();
    if(!(R_pivots.minCoeff()
         > Scalar{static_cast<int>(NU)} * std::numeric_limits<Scalar>::epsilon() * R_pivots.maxCoeff()))
        return ctrlpp::unexpected(care_error::singular_r);

    const MatU R_inv = R_ldlt.solve(MatU::Identity());

    Mat2N H;
    H.template block<nx, nx>(0, 0) = A;
    H.template block<nx, nx>(0, nx) = -B * R_inv * B.transpose();
    H.template block<nx, nx>(nx, 0) = -Q;
    H.template block<nx, nx>(nx, nx) = -A.transpose();

    // Still reachable after the rank test above, and for a different reason than
    // that test covers: an R that passes it but sits close to the threshold gives a
    // large R^{-1}, and B R^{-1} B^T can overflow for a large enough B. That is an
    // assembled quantity leaving the representable range, which is what the
    // enumerator names -- not a statement about the arguments the caller passed.
    if(!H.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    auto result = detail::care_solve_from_hamiltonian<Scalar, NX, Method>(H);
    if(!result)
        return ctrlpp::unexpected(result.error());

    return Eigen::Matrix<Scalar, int(NU), int(NX)>{(R_inv * B.transpose() * result->P).eval()};
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
                        const Eigen::Matrix<Scalar, int(NU), int(NU)>& R) -> ctrlpp::expected<lqi_result<Scalar, NX, NU, NY>, dare_error>
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int ny = static_cast<int>(NY);
    constexpr int nu = static_cast<int>(NU);

    auto K_aug_opt = lqr_gain<Scalar, NX + NY, NU>(A_aug, B_aug, Q_aug, R);
    if(!K_aug_opt)
        return ctrlpp::unexpected(K_aug_opt.error());

    auto& K_aug = *K_aug_opt;
    lqi_result<Scalar, NX, NU, NY> result;
    result.Kx = K_aug.template block<nu, nx>(0, 0);
    result.Ki = K_aug.template block<nu, ny>(0, nx);
    return result;
}

}

/// LQI gain: augments state with integral of tracking error.
///
/// Augmented system: A_aug = [[A, 0], [-C, I]], B_aug = [[B], [0]].
/// Returns lqi_result with partitioned Kx (NU x NX) and Ki (NU x NY), or the
/// augmented Riccati solve's `dare_error` verbatim.
///
/// @cite anderson1990 -- Anderson & Moore, "Optimal Control: Linear Quadratic Methods", 1990, Ch. 9 (integral action)
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
auto lqi_gain(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
              const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
              const Eigen::Matrix<Scalar, int(NY), int(NX)>& C,
              const Eigen::Matrix<Scalar, int(NX + NY), int(NX + NY)>& Q_aug,
              const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
    -> ctrlpp::expected<lqi_result<Scalar, NX, NU, NY>, dare_error>
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
template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU>
class lqr
{
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
template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU>
class lqr_time_varying
{
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
