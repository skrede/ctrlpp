#ifndef HPP_GUARD_CTRLPP_CONTROL_LQR_H
#define HPP_GUARD_CTRLPP_CONTROL_LQR_H

// @cite anderson1990 -- Anderson & Moore, "Optimal Control: Linear Quadratic Methods", 1990

#include "ctrlpp/control/care.h"
#include "ctrlpp/control/dare.h"
#include "ctrlpp/control/lqr_types.h"

#include "ctrlpp/control/detail/lqr_discrete.h"
#include "ctrlpp/control/detail/lqr_continuous.h"
#include "ctrlpp/control/detail/care_hamiltonian.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include <Eigen/Dense>

#include <span>
#include <vector>
#include <cstddef>
#include <utility>

namespace ctrlpp
{

// Returns the solve's own verified gain: re-forming it from P would take a second, unthresholded factorization of
// R + B'PB at the caller's scale, which can hand back a least-squares gain or overflow on a certified pose.
template <typename Scalar, std::size_t NX, std::size_t NU>
ctrlpp::expected<Eigen::Matrix<Scalar, int(NU), int(NX)>, dare_error> lqr_gain(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
{
    auto P_result = dare<Scalar, NX, NU>(A, B, Q, R);
    if(!P_result)
        return ctrlpp::unexpected(P_result.error());
    return P_result->K;
}

template <typename Scalar, std::size_t NX, std::size_t NU>
ctrlpp::expected<Eigen::Matrix<Scalar, int(NU), int(NX)>, dare_error> lqr_gain(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R, const Eigen::Matrix<Scalar, int(NX), int(NU)>& N)
{
    auto P_result = dare<Scalar, NX, NU>(A, B, Q, R, N);
    if(!P_result)
        return ctrlpp::unexpected(P_result.error());
    return P_result->K;
}

// Forms R^{-1} once and builds the Hamiltonian itself, so the one factorization of R serves both the solve and the gain.
template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, detail::care_solve_method Method = detail::sign_function_care_method>
ctrlpp::expected<Eigen::Matrix<Scalar, int(NU), int(NX)>, care_error> lqr_gain_continuous(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
        const Eigen::Matrix<Scalar, int(NX), int(NU)>& B, const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R, Method /*method_tag*/ = {})
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

    if(!A.allFinite() || !B.allFinite() || !Q.allFinite() || !R.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    // The same equilibration as care(), so the two entry points answer the same poses. Nothing is unscaled
    // afterwards: the gain is homogeneous of degree zero in the common weight scale.
    auto Q_scaled = Q;
    auto R_scaled = R;
    if(!detail::equilibrate_care_weights<Scalar, NX, NU>(detail::care_weight_scale<Scalar, NX, NU>(Q, R), Q_scaled, R_scaled))
        return ctrlpp::unexpected(care_error::non_finite_input);

    Eigen::Matrix<Scalar, int(NU), int(NU)> R_inv;
    const auto H = detail::lqr_continuous_hamiltonian<Scalar, NX, NU>(A, B, Q_scaled, R_scaled, R_inv);
    if(!H)
        return ctrlpp::unexpected(H.error());
    auto result = detail::care_solve_from_hamiltonian<Scalar, NX, Method>(*H);
    if(!result)
        return ctrlpp::unexpected(result.error());

    return Eigen::Matrix<Scalar, int(NU), int(NX)>{(R_inv * B.transpose() * result->P).eval()};
}

template <typename Scalar, std::size_t NX, std::size_t NU>
std::vector<Eigen::Matrix<Scalar, int(NU), int(NX)>> lqr_finite(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R, const Eigen::Matrix<Scalar, int(NX), int(NX)>& Qf, std::size_t horizon)
{
    std::vector<Eigen::Matrix<Scalar, int(NU), int(NX)>> gains(horizon);
    Eigen::Matrix<Scalar, int(NX), int(NX)> P = Qf;
    for(std::size_t i = horizon; i > 0; --i)
        gains[i - 1] = detail::lqr_riccati_step<Scalar, NX, NU>(A, B, Q, R, P);
    return gains;
}

template <typename Scalar, std::size_t NX, std::size_t NU>
ctrlpp::expected<std::vector<Eigen::Matrix<Scalar, int(NU), int(NX)>>, lqr_error> lqr_tv_gains(const std::vector<Eigen::Matrix<Scalar, int(NX), int(NX)>>& As,
        const std::vector<Eigen::Matrix<Scalar, int(NX), int(NU)>>& Bs, const std::vector<Eigen::Matrix<Scalar, int(NX), int(NX)>>& Qs,
        const std::vector<Eigen::Matrix<Scalar, int(NU), int(NU)>>& Rs, const Eigen::Matrix<Scalar, int(NX), int(NX)>& Qf)
{
    if(auto checked = detail::check_lqr_tv_sequences<Scalar, NX, NU>(As, Bs, Qs, Rs, Qf); !checked)
        return ctrlpp::unexpected(checked.error());

    std::vector<Eigen::Matrix<Scalar, int(NU), int(NX)>> gains(As.size());
    Eigen::Matrix<Scalar, int(NX), int(NX)> P = Qf;
    for(std::size_t i = As.size(); i > 0; --i)
        gains[i - 1] = detail::lqr_riccati_step<Scalar, NX, NU>(As[i - 1], Bs[i - 1], Qs[i - 1], Rs[i - 1], P);
    return gains;
}

// The state is augmented with the integral of the tracking error: A_aug = [[A, 0], [-C, I]], B_aug = [[B], [0]].
// @cite anderson1990 -- Anderson & Moore, "Optimal Control: Linear Quadratic Methods", 1990, Ch. 9 (integral action)
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
ctrlpp::expected<lqi_result<Scalar, NX, NU, NY>, dare_error> lqi_gain(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
        const Eigen::Matrix<Scalar, int(NY), int(NX)>& C, const Eigen::Matrix<Scalar, int(NX + NY), int(NX + NY)>& Q_aug, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
{
    auto [A_aug, B_aug] = detail::build_lqi_augmented_system<Scalar, NX, NU, NY>(A, B, C);
    auto K_aug = lqr_gain<Scalar, NX + NY, NU>(A_aug, B_aug, Q_aug, R);
    if(!K_aug)
        return ctrlpp::unexpected(K_aug.error());
    return detail::partition_lqi_gain<Scalar, NX, NU, NY>(*K_aug);
}

template <typename Scalar, std::size_t NX, std::size_t NU>
ctrlpp::expected<Scalar, lqr_error> lqr_cost(std::span<const Eigen::Matrix<Scalar, int(NX), 1>> xs, std::span<const Eigen::Matrix<Scalar, int(NU), 1>> us,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
{
    if(xs.size() != us.size() && xs.size() != us.size() + 1)
        return ctrlpp::unexpected(lqr_error::trajectory_length_mismatch);
    if(!detail::lqr_cost_operands_finite<Scalar, NX, NU>(xs, us, Q, R))
        return ctrlpp::unexpected(lqr_error::non_finite_input);
    return detail::lqr_trajectory_cost<Scalar, NX, NU>(xs, us, Q, R);
}

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU>
class lqr
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

public:
    using gain_type = Eigen::Matrix<Scalar, int(NU), int(NX)>;
    using state_type = Eigen::Matrix<Scalar, int(NX), 1>;
    using input_type = Eigen::Matrix<Scalar, int(NU), 1>;

    explicit lqr(gain_type K)
        : m_K(std::move(K))
    {
    }

    input_type compute(const state_type& x) const { return (-m_K * x).eval(); }

    const gain_type& gain() const { return m_K; }

private:
    gain_type m_K;
};

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU>
class lqr_time_varying
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

public:
    using gain_type = Eigen::Matrix<Scalar, int(NU), int(NX)>;
    using state_type = Eigen::Matrix<Scalar, int(NX), 1>;
    using input_type = Eigen::Matrix<Scalar, int(NU), 1>;

    explicit lqr_time_varying(std::vector<gain_type> gains)
        : m_gains(std::move(gains))
    {
    }

    ctrlpp::expected<input_type, lqr_error> compute(const state_type& x, std::size_t k) const
    {
        if(k >= m_gains.size())
            return ctrlpp::unexpected(lqr_error::step_out_of_range);
        return (-m_gains[k] * x).eval();
    }

    ctrlpp::expected<gain_type, lqr_error> gain(std::size_t k) const
    {
        if(k >= m_gains.size())
            return ctrlpp::unexpected(lqr_error::step_out_of_range);
        return m_gains[k];
    }

    std::size_t horizon() const { return m_gains.size(); }

private:
    std::vector<gain_type> m_gains;
};

}

#endif
