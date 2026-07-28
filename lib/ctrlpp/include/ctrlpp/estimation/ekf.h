#ifndef HPP_GUARD_CTRLPP_ESTIMATION_EKF_H
#define HPP_GUARD_CTRLPP_ESTIMATION_EKF_H

/// @brief Extended Kalman Filter with analytical/numerical Jacobian dispatch.
///
/// @cite simon2006 -- Simon, "Optimal State Estimation", 2006, Ch. 13

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/model/measurement_model.h"
#include "ctrlpp/model/differentiable_dynamics.h"
#include "ctrlpp/model/differentiable_measurement.h"

#include "ctrlpp/detail/covariance_ops.h"
#include "ctrlpp/detail/numerical_diff.h"

#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/estimation_types.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <utility>

namespace ctrlpp
{

/// @brief Structured failure modes of an `ekf` measurement update.
///
/// Each enumerator is an exact domain condition, not a tuning preference: a
/// non-finite operand makes every downstream product non-finite, so the step
/// cannot produce an estimate at all.
///
///  * non_finite_state       : the carried state estimate is already non-finite
///                             when the step begins. It is reported ahead of the
///                             measurement because the fault is upstream of it,
///                             and because the measurement Jacobian is evaluated
///                             AT the carried state, so a poisoned state makes
///                             the linearization meaningless before the
///                             measurement is ever used.
///  * non_finite_covariance  : the carried covariance is already non-finite when
///                             the step begins. A distinguishable cause from a
///                             non-finite state because the covariance recursion
///                             is driven by the linearized dynamics and by Q and
///                             R, never by the measurement.
///  * non_finite_measurement : the supplied measurement vector has a non-finite
///                             component. The gain carries it into the state,
///                             which is the filter's carried memory, so a single
///                             such sample destroys the estimate permanently.
enum class ekf_update_error
{
    non_finite_state,
    non_finite_covariance,
    non_finite_measurement,
};

/// @brief Persistent state-health status of an `ekf`.
///
/// A per-call result cannot answer whether the carried estimate is still
/// degraded from a step several samples ago, because that question outlives the
/// call. The status latches: it never returns to a lower rung on its own.
///
///  * ok                  : every step so far began from a finite estimate.
///  * non_finite_estimate : a step found the carried state or the carried
///                          covariance already non-finite. `predict` does not
///                          reject its input, so this is how a poisoned control
///                          vector or a dynamics model that returned a
///                          non-finite state becomes visible. A rejected
///                          measurement does NOT set it: the rejection mutates
///                          nothing, so it leaves the filter healthy.
enum class ekf_health
{
    ok,
    non_finite_estimate,
};

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct ekf_config
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NY > 0, "Output dimension NY must be positive");
    Matrix<Scalar, NX, NX> Q{Matrix<Scalar, NX, NX>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Vector<Scalar, NX> x0{Vector<Scalar, NX>::Zero()};
    Matrix<Scalar, NX, NX> P0{Matrix<Scalar, NX, NX>::Identity()};
    Scalar numerical_eps{std::cbrt(std::numeric_limits<Scalar>::epsilon())};
};

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY, typename Dynamics, typename Measurement>
    requires dynamics_model<Dynamics, Scalar, NX, NU> && measurement_model<Measurement, Scalar, NX, NY>
class ekf
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NY > 0, "Output dimension NY must be positive");

    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);
    static constexpr int ny = static_cast<int>(NY);

public:
    using observer_tag = struct ekf_tag;
    using state_vector_t = Vector<Scalar, NX>;
    using input_vector_t = Vector<Scalar, NU>;
    using output_vector_t = Vector<Scalar, NY>;
    using cov_matrix_t = Matrix<Scalar, NX, NX>;
    using meas_cov_matrix_t = Matrix<Scalar, NY, NY>;

    /// @brief Fallible factory, and the only way to originate a filter.
    ///
    /// Rejects a configuration whose process noise, measurement noise, initial
    /// state or initial covariance carries a non-finite entry, naming which of
    /// the four it is. Such a configuration is a mistake the caller made before
    /// the filter ever ran: an infinite Q makes the predicted covariance
    /// infinite, the gain a ratio of infinities and the estimate non-finite at
    /// the first step, which surfaces the fault as far as possible from where it
    /// was introduced.
    ///
    /// Finiteness is the domain condition and the whole of it. A covariance that
    /// is merely ill-conditioned -- entries many orders of magnitude apart, or a
    /// singular P0 -- is a legitimately posed problem and is accepted. Rejecting
    /// it would convert a numerical-behavior question into a domain violation
    /// and refuse configurations that work.
    ///
    /// `numerical_eps` is NOT validated here. It is the finite-difference step
    /// used only when the model is not analytically differentiable, and it
    /// carries its own domain condition (finite and strictly positive, since the
    /// central-difference stencil divides by it); converting that condition is
    /// separately owned work.
    static auto create(Dynamics dynamics, Measurement measurement, ekf_config<Scalar, NX, NU, NY> config) -> ctrlpp::expected<ekf, filter_error>
    {
        if(const auto valid = detail::validate_filter_configuration(config.Q, config.R, config.x0, config.P0); !valid)
            return ctrlpp::unexpected(valid.error());
        return ekf{validated_tag{}, std::move(dynamics), std::move(measurement), std::move(config)};
    }

    void predict(const input_vector_t& u)
    {
        auto x_prev = m_x;
        propagate_state(x_prev, u);
        propagate_covariance(x_prev, u);
    }

    /// @brief Update state and covariance with a measurement.
    ///
    /// The step is rejected before any member is assigned when the carried
    /// estimate or the measurement is non-finite, so a rejected step leaves the
    /// state, the covariance, the innovation and the NIS bitwise unchanged and
    /// the caller may retry with the next sample.
    ///
    /// `predict` is deliberately not fallible: its input is a control vector the
    /// caller already commanded and owns, and rejecting it would leave the
    /// filter with no propagation for a step the plant did take. A prediction
    /// that poisons the state is instead reported by `health()`, which the next
    /// update latches.
    auto update(const output_vector_t& z) -> ctrlpp::expected<void, ekf_update_error>
    {
        if(const auto step = check_step(z); !step)
            return ctrlpp::unexpected(latch_health(step.error()));

        auto z_pred = m_measurement(m_x);
        m_innovation = (z - z_pred).eval();

        auto H = compute_measurement_jacobian();
        auto S = compute_innovation_covariance(H);
        auto K = compute_kalman_gain(H, S);

        apply_state_correction(K, z, z_pred);
        update_covariance(K, H);

        m_nis = (m_innovation.transpose() * S.colPivHouseholderQr().solve(m_innovation))(0, 0);
        return {};
    }

    const state_vector_t& state() const { return m_x; }

    const cov_matrix_t& covariance() const { return m_P; }

    const output_vector_t& innovation() const { return m_innovation; }

    /// @brief Normalized Innovation Squared: innovation^T S^{-1} innovation
    /// (chi-square distributed with dof = NY under a consistent filter).
    Scalar nis() const { return m_nis; }

    /// @brief Report whether the carried estimate is still degraded from an
    /// earlier step. Latches; a rejected measurement does not set it.
    ekf_health health() const { return m_health; }

private:
    /// @brief Tag selecting the non-validating constructor reserved for
    /// `create`, which is what makes the factory the only public path and the
    /// validation impossible to bypass.
    struct validated_tag
    {
    };

    ekf(validated_tag, Dynamics dynamics, Measurement measurement, ekf_config<Scalar, NX, NU, NY> config)
        : m_eps{config.numerical_eps}
        , m_dynamics{std::move(dynamics)}
        , m_P{std::move(config.P0)}
        , m_Q{std::move(config.Q)}
        , m_x{std::move(config.x0)}
        , m_R{std::move(config.R)}
        , m_measurement{std::move(measurement)}
        , m_innovation{output_vector_t::Zero()}
    {
    }

    /// @brief Classify a step's operands without touching a single member.
    ///
    /// The order is the severity order documented on `ekf_update_error`: the
    /// carried estimate first, the supplied measurement last. The cost is one
    /// finiteness scan of each operand -- NX + NX*NX + NY reads, no branches on
    /// data and no allocation, all dimensions being compile-time constants.
    auto check_step(const output_vector_t& z) const -> ctrlpp::expected<void, ekf_update_error>
    {
        if(!m_x.allFinite())
            return ctrlpp::unexpected(ekf_update_error::non_finite_state);
        if(!m_P.allFinite())
            return ctrlpp::unexpected(ekf_update_error::non_finite_covariance);
        if(!z.allFinite())
            return ctrlpp::unexpected(ekf_update_error::non_finite_measurement);
        return {};
    }

    /// @brief Latch the persistent status for the faults that describe the
    /// carried estimate, and pass the fault through so the caller still receives
    /// the specific cause on the failure channel.
    auto latch_health(ekf_update_error fault) -> ekf_update_error
    {
        if(fault != ekf_update_error::non_finite_measurement)
            m_health = ekf_health::non_finite_estimate;
        return fault;
    }

    /// @brief Propagate state through dynamics model.
    void propagate_state(const state_vector_t& x_prev, const input_vector_t& u)
    {
        m_x = m_dynamics(x_prev, u);
    }

    /// @brief Propagate covariance through linearized dynamics: P = F*P*F^T + Q.
    ///
    /// @cite simon2006 -- Simon, "Optimal State Estimation", 2006, Ch. 13, Eq. 13.1
    void propagate_covariance(const state_vector_t& x_prev, const input_vector_t& u)
    {
        Matrix<Scalar, NX, NX> F;
        if constexpr(differentiable_dynamics<Dynamics, Scalar, NX, NU>)
            F = m_dynamics.jacobian_x(x_prev, u);
        else
            F = detail::numerical_jacobian_x<Scalar, NX, NU>(m_dynamics, x_prev, u, m_eps);

        m_P = detail::symmetrize((F * m_P * F.transpose() + m_Q).eval());
    }

    /// @brief Compute measurement Jacobian H at current state (analytical or numerical).
    auto compute_measurement_jacobian() const -> Matrix<Scalar, NY, NX>
    {
        if constexpr(differentiable_measurement<Measurement, Scalar, NX, NY>)
            return m_measurement.jacobian(m_x);
        else
            return detail::numerical_jacobian_h<Scalar, NX, NY>(m_measurement, m_x, m_eps);
    }

    /// @brief Compute innovation covariance: S = H*P*H^T + R.
    ///
    /// @cite simon2006 -- Simon, "Optimal State Estimation", 2006, Ch. 13, Eq. 13.2
    auto compute_innovation_covariance(const Matrix<Scalar, NY, NX>& H) const -> meas_cov_matrix_t
    {
        return (H * m_P * H.transpose() + m_R).eval();
    }

    /// @brief Compute Kalman gain: K = P*H^T*S^{-1} via column-pivoting QR solve.
    ///
    /// @cite simon2006 -- Simon, "Optimal State Estimation", 2006, Ch. 13, Eq. 13.2
    auto compute_kalman_gain(const Matrix<Scalar, NY, NX>& H, const meas_cov_matrix_t& S) const -> Eigen::Matrix<Scalar, nx, ny>
    {
        Eigen::Matrix<Scalar, ny, nx> KT_solved = S.transpose().colPivHouseholderQr().solve(H * m_P);
        return KT_solved.transpose().eval();
    }

    /// @brief Apply state correction: x += K * innovation.
    void apply_state_correction(const Eigen::Matrix<Scalar, nx, ny>& K, const output_vector_t& /*z*/, const output_vector_t& /*z_pred*/)
    {
        m_x = (m_x + K * m_innovation).eval();
    }

    /// @brief Joseph-form covariance update: P = (I-KH)*P*(I-KH)^T + K*R*K^T.
    ///
    /// @cite simon2006 -- Simon, "Optimal State Estimation", 2006, Ch. 13, Eq. 13.3
    void update_covariance(const Eigen::Matrix<Scalar, nx, ny>& K, const Matrix<Scalar, NY, NX>& H)
    {
        cov_matrix_t IKH = cov_matrix_t::Identity() - K * H;
        m_P = detail::symmetrize((IKH * m_P * IKH.transpose() + K * m_R * K.transpose()).eval());
    }

    Scalar m_eps;
    Scalar m_nis{0};
    Dynamics m_dynamics;
    cov_matrix_t m_P;
    cov_matrix_t m_Q;
    state_vector_t m_x;
    meas_cov_matrix_t m_R;
    Measurement m_measurement;
    output_vector_t m_innovation;
    ekf_health m_health{ekf_health::ok};
};

namespace detail
{

struct ekf_sa_dynamics
{
    Vector<double, 2> operator()(const Vector<double, 2>&, const Vector<double, 1>&) const { return Vector<double, 2>::Zero(); }
};

struct ekf_sa_measurement
{
    Vector<double, 1> operator()(const Vector<double, 2>&) const { return Vector<double, 1>::Zero(); }
};

}

static_assert(ObserverPolicy<ekf<double, 2, 1, 1, detail::ekf_sa_dynamics, detail::ekf_sa_measurement>>);
static_assert(CovarianceObserver<ekf<double, 2, 1, 1, detail::ekf_sa_dynamics, detail::ekf_sa_measurement>>);

}

#endif
