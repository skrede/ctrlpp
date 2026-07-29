#ifndef HPP_GUARD_CTRLPP_ESTIMATION_KALMAN_H
#define HPP_GUARD_CTRLPP_ESTIMATION_KALMAN_H

/// @brief Linear discrete-time Kalman filter with Joseph-form covariance update.
///
/// @cite kalman1960 -- Kalman, "A New Approach to Linear Filtering and Prediction Problems", 1960

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/model/state_space.h"

#include "ctrlpp/detail/covariance_ops.h"

#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/estimation_types.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <utility>

namespace ctrlpp
{

/// @brief Structured failure modes of a `kalman_filter` measurement update.
///
/// Each enumerator is an exact domain condition, not a tuning preference: a
/// non-finite operand makes every downstream product non-finite, so the step
/// cannot produce an estimate at all.
///
///  * non_finite_state       : the carried state estimate is already non-finite
///                             when the step begins. The fault is upstream of
///                             the measurement, so it is reported ahead of it --
///                             a caller told "your measurement is bad" would
///                             replace a working sensor while the real fault
///                             sits in the prediction that poisoned the state.
///  * non_finite_covariance  : the carried covariance is already non-finite when
///                             the step begins. A distinguishable cause from a
///                             non-finite state because the covariance recursion
///                             is driven by the model and by Q and R, never by
///                             the measurement, so the caller fixes a different
///                             input.
///  * non_finite_measurement : the supplied measurement vector has a non-finite
///                             component. The gain carries it into the state,
///                             which is the filter's carried memory, so a single
///                             such sample destroys the estimate permanently.
enum class kalman_update_error
{
    non_finite_state,
    non_finite_covariance,
    non_finite_measurement,
    non_finite_result,
};

/// @brief Persistent state-health status of a `kalman_filter`.
///
/// A per-call result cannot answer whether the carried estimate is still
/// degraded from a step several samples ago, because that question outlives the
/// call. The status latches: it never returns to a lower rung on its own.
///
///  * ok                  : every step so far began from a finite estimate.
///  * non_finite_estimate : a step found the carried state or the carried
///                          covariance already non-finite. `predict` does not
///                          reject its input, so this is how a poisoned control
///                          vector or a non-finite model becomes visible. A
///                          rejected measurement does NOT set it: the rejection
///                          mutates nothing, so it leaves the filter healthy.
enum class kalman_health
{
    ok,
    non_finite_estimate,
};

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct kalman_config
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NY > 0, "Output dimension NY must be positive");
    Matrix<Scalar, NX, NX> Q{Matrix<Scalar, NX, NX>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Vector<Scalar, NX> x0{Vector<Scalar, NX>::Zero()};
    Matrix<Scalar, NX, NX> P0{Matrix<Scalar, NX, NX>::Identity()};
};

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
class kalman_filter
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NY > 0, "Output dimension NY must be positive");

    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);
    static constexpr int ny = static_cast<int>(NY);

public:
    using observer_tag = struct kalman_tag;
    using state_vector_t = Eigen::Matrix<Scalar, nx, 1>;
    using input_vector_t = Eigen::Matrix<Scalar, nu, 1>;
    using output_vector_t = Eigen::Matrix<Scalar, ny, 1>;
    using cov_matrix_t = Eigen::Matrix<Scalar, nx, nx>;
    using meas_cov_matrix_t = Eigen::Matrix<Scalar, ny, ny>;
    using system_t = discrete_state_space<Scalar, NX, NU, NY>;

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
    static auto create(system_t sys, kalman_config<Scalar, NX, NU, NY> config) -> ctrlpp::expected<kalman_filter, filter_error>
    {
        if(const auto valid = detail::validate_filter_configuration(config.Q, config.R, config.x0, config.P0); !valid)
            return ctrlpp::unexpected(valid.error());
        return kalman_filter{validated_tag{}, std::move(sys), std::move(config)};
    }

    /// @brief Predict state and covariance one step forward.
    void predict(const input_vector_t& u)
    {
        m_P_post_prev = m_P;
        propagate_state(u);
        propagate_covariance();
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
    auto update(const output_vector_t& z) -> ctrlpp::expected<void, kalman_update_error>
    {
        if(const auto step = check_step(z); !step)
            return ctrlpp::unexpected(latch_health(step.error()));

        auto const previous_x = m_x;
        auto const previous_P = m_P;
        auto const previous_innovation = m_innovation;
        auto const previous_nis = m_nis_value;
        compute_innovation(z);
        auto S = compute_innovation_covariance();
        auto K = compute_kalman_gain(S);

        apply_state_correction(K);
        update_covariance(K);
        compute_nis(S);
        if(!m_x.allFinite() || !m_P.allFinite()
            || !m_innovation.allFinite() || !std::isfinite(m_nis_value))
        {
            m_x = previous_x;
            m_P = previous_P;
            m_innovation = previous_innovation;
            m_nis_value = previous_nis;
            return ctrlpp::unexpected(kalman_update_error::non_finite_result);
        }
        return {};
    }

    const state_vector_t& state() const { return m_x; }

    const cov_matrix_t& covariance() const { return m_P; }

    const output_vector_t& innovation() const { return m_innovation; }

    /// @brief Normalized Innovation Squared: innovation^T S^{-1} innovation
    /// (chi-square distributed with dof = NY under a consistent filter).
    Scalar nis() const { return m_nis_value; }

    /// @brief Report whether the carried estimate is still degraded from an
    /// earlier step. Latches; a rejected measurement does not set it.
    kalman_health health() const { return m_health; }

    bool is_steady_state(Scalar tol = Scalar{1e-10}) const
    {
        Scalar p_norm = m_P.norm();
        if(p_norm < std::numeric_limits<Scalar>::epsilon())
            return true;
        return (m_P - m_P_post_prev).norm() / p_norm < tol;
    }

    void reset_covariance(const cov_matrix_t& P0)
    {
        m_P = P0;
        m_P_post_prev = P0;
    }

    void set_model(system_t sys) { m_sys = std::move(sys); }

    void set_noise(cov_matrix_t Q, meas_cov_matrix_t R)
    {
        m_Q = std::move(Q);
        m_R = std::move(R);
    }

private:
    /// @brief Tag selecting the non-validating constructor reserved for
    /// `create`, which is what makes the factory the only public path and the
    /// validation impossible to bypass.
    struct validated_tag
    {
    };

    kalman_filter(validated_tag, system_t sys, kalman_config<Scalar, NX, NU, NY> config)
        : m_sys{std::move(sys)}, m_Q{std::move(config.Q)}, m_P{std::move(config.P0)}, m_P_post_prev{m_P}, m_x{std::move(config.x0)}, m_R{std::move(config.R)}, m_innovation{output_vector_t::Zero()}
    {
    }

    /// @brief Classify a step's operands without touching a single member.
    ///
    /// The order is the severity order documented on `kalman_update_error`: the
    /// carried estimate first, the supplied measurement last. The cost is one
    /// finiteness scan of each operand -- NX + NX*NX + NY reads, no branches on
    /// data and no allocation, all dimensions being compile-time constants.
    auto check_step(const output_vector_t& z) const -> ctrlpp::expected<void, kalman_update_error>
    {
        if(!m_x.allFinite())
            return ctrlpp::unexpected(kalman_update_error::non_finite_state);
        if(!m_P.allFinite())
            return ctrlpp::unexpected(kalman_update_error::non_finite_covariance);
        if(!z.allFinite())
            return ctrlpp::unexpected(kalman_update_error::non_finite_measurement);
        return {};
    }

    /// @brief Latch the persistent status for the faults that describe the
    /// carried estimate, and pass the fault through so the caller still receives
    /// the specific cause on the failure channel.
    auto latch_health(kalman_update_error fault) -> kalman_update_error
    {
        if(fault != kalman_update_error::non_finite_measurement)
            m_health = kalman_health::non_finite_estimate;
        return fault;
    }

    /// @brief Propagate state: x = A*x + B*u.
    ///
    /// @cite kalman1960 -- Kalman, "A New Approach to Linear Filtering and Prediction Problems", 1960
    void propagate_state(const input_vector_t& u)
    {
        m_x = (m_sys.A * m_x + m_sys.B * u).eval();
    }

    /// @brief Propagate covariance: P = A*P*A^T + Q.
    ///
    /// @cite kalman1960 -- Kalman, "A New Approach to Linear Filtering and Prediction Problems", 1960
    void propagate_covariance()
    {
        m_P = (m_sys.A * m_P * m_sys.A.transpose() + m_Q).eval();
    }

    /// @brief Compute innovation: y = z - C*x.
    void compute_innovation(const output_vector_t& z)
    {
        m_innovation = (z - m_sys.C * m_x).eval();
    }

    /// @brief Compute innovation covariance: S = C*P*C^T + R.
    auto compute_innovation_covariance() const -> meas_cov_matrix_t
    {
        return (m_sys.C * m_P * m_sys.C.transpose() + m_R).eval();
    }

    /// @brief Compute Kalman gain: K = P*C^T*S^{-1} via column-pivoting QR solve.
    ///
    /// @cite kalman1960 -- Kalman, "A New Approach to Linear Filtering and Prediction Problems", 1960
    auto compute_kalman_gain(const meas_cov_matrix_t& S) const -> Eigen::Matrix<Scalar, nx, ny>
    {
        Eigen::Matrix<Scalar, ny, nx> CP = m_sys.C * m_P;
        Eigen::Matrix<Scalar, ny, nx> KT_solved = S.transpose().colPivHouseholderQr().solve(CP);
        return KT_solved.transpose().eval();
    }

    /// @brief Apply state correction: x += K * innovation.
    void apply_state_correction(const Eigen::Matrix<Scalar, nx, ny>& K)
    {
        m_x = (m_x + K * m_innovation).eval();
    }

    /// @brief Joseph-form covariance update: P = (I-KC)*P*(I-KC)^T + K*R*K^T, symmetrized.
    ///
    /// @cite kalman1960 -- Kalman, "A New Approach to Linear Filtering and Prediction Problems", 1960
    void update_covariance(const Eigen::Matrix<Scalar, nx, ny>& K)
    {
        cov_matrix_t IKC = cov_matrix_t::Identity() - K * m_sys.C;
        m_P = detail::symmetrize((IKC * m_P * IKC.transpose() + K * m_R * K.transpose()).eval());
    }

    /// @brief Compute Normalized Innovation Squared: innovation^T * S^{-1} * innovation
    /// (chi-square distributed with dof = NY).
    void compute_nis(const meas_cov_matrix_t& S)
    {
        output_vector_t Sinv_z = S.colPivHouseholderQr().solve(m_innovation).eval();
        m_nis_value = (m_innovation.transpose() * Sinv_z)(0, 0);
    }

    Scalar m_nis_value{0};
    system_t m_sys;
    cov_matrix_t m_Q;
    cov_matrix_t m_P;
    cov_matrix_t m_P_post_prev;
    state_vector_t m_x;
    meas_cov_matrix_t m_R;
    output_vector_t m_innovation;
    kalman_health m_health{kalman_health::ok};
};

static_assert(ObserverPolicy<kalman_filter<double, 2, 1, 1>>);
static_assert(CovarianceObserver<kalman_filter<double, 2, 1, 1>>);

}

#endif
