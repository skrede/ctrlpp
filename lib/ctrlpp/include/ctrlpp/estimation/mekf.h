#ifndef HPP_GUARD_CTRLPP_ESTIMATION_MEKF_H
#define HPP_GUARD_CTRLPP_ESTIMATION_MEKF_H

/// @brief Multiplicative Extended Kalman Filter for attitude estimation.
///
/// Two-track state: 7D nominal (quaternion + bias) with 6D tangent-space
/// error-state covariance. The mandatory post-update covariance reset via
/// frame-change Jacobian G is the key correctness concern.
///
/// @cite markley2003 -- Markley, "Attitude Error Representations for Kalman Filtering", 2003
/// @cite crassidis2003 -- Crassidis & Markley, "Unscented Filtering for Spacecraft Attitude Estimation", J. Guidance Control Dyn 26(4), 2003

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/lie/so3.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/detail/covariance_ops.h"
#include "ctrlpp/detail/numerical_mekf_diff.h"

#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/estimation_types.h"

#include <Eigen/Geometry>

#include <cmath>
#include <limits>
#include <cstddef>
#include <utility>

namespace ctrlpp
{

template <typename M, typename Scalar, std::size_t NB, std::size_t NY>
concept mekf_measurement_model = requires(const M& m, const Eigen::Quaternion<Scalar>& q, const Vector<Scalar, NB>& b) {
    { m(q, b) } -> std::convertible_to<Vector<Scalar, NY>>;
};

template <typename M, typename Scalar, std::size_t NB, std::size_t NY>
concept differentiable_mekf_measurement = mekf_measurement_model<M, Scalar, NB, NY> && requires(const M& m, const Eigen::Quaternion<Scalar>& q, const Vector<Scalar, NB>& b) {
    { m.jacobian(q, b) } -> std::convertible_to<Matrix<Scalar, NY, 3 + NB>>;
};

/// @brief Structured failure modes of a `mekf` measurement update.
///
/// Each enumerator is an exact domain condition, not a tuning preference: a
/// non-finite operand makes every downstream product non-finite, so the step
/// cannot produce an estimate at all.
///
///  * non_finite_state       : the carried nominal state -- the attitude
///                             quaternion or the bias -- is already non-finite
///                             when the step begins. It is reported ahead of the
///                             measurement because the fault is upstream of it,
///                             and because the measurement and its Jacobian are
///                             evaluated AT that nominal state.
///  * non_finite_covariance  : the carried error-state covariance is already
///                             non-finite when the step begins. A distinguishable
///                             cause from a non-finite nominal state because the
///                             covariance recursion is driven by the propagation
///                             Jacobian and by Q and R, never by the measurement.
///  * non_finite_measurement : the supplied measurement vector has a non-finite
///                             component. The gain carries it into the
///                             multiplicative correction, so a single such
///                             sample makes the attitude quaternion non-finite
///                             and no later normalization recovers it.
enum class mekf_update_error
{
    non_finite_state,
    non_finite_covariance,
    non_finite_measurement,
};

/// @brief Persistent state-health status of a `mekf`.
///
/// A per-call result cannot answer whether the carried estimate is still
/// degraded from a step several samples ago, because that question outlives the
/// call. The status latches: it never returns to a lower rung on its own.
///
///  * ok                  : every step so far began from a finite estimate.
///  * non_finite_estimate : a step found the carried nominal state or the
///                          carried covariance already non-finite. `predict`
///                          does not reject its input, so this is how a poisoned
///                          gyro rate or a non-finite timestep becomes visible.
///                          A rejected measurement does NOT set it: the
///                          rejection mutates nothing, so it leaves the filter
///                          healthy.
enum class mekf_health
{
    ok,
    non_finite_estimate,
};

template <ctrlpp_floating_scalar Scalar, std::size_t NB, std::size_t NY>
struct mekf_config
{
    static_assert(NB >= 3, "Bias dimension NB must be at least 3: the propagation subtracts the leading three bias elements from the gyro rate");
    static_assert(NY > 0, "Output dimension NY must be positive");
    static constexpr std::size_t NE = 3 + NB;
    Matrix<Scalar, NE, NE> Q{Matrix<Scalar, NE, NE>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Eigen::Quaternion<Scalar> q0{Eigen::Quaternion<Scalar>::Identity()};
    Vector<Scalar, NB> b0{Vector<Scalar, NB>::Zero()};
    Matrix<Scalar, NE, NE> P0{Matrix<Scalar, NE, NE>::Identity()};
    Scalar dt{Scalar{0.01}};
    Scalar numerical_eps{std::cbrt(std::numeric_limits<Scalar>::epsilon())};
};

template <ctrlpp_floating_scalar Scalar, std::size_t NB, std::size_t NY, typename Measurement>
    requires mekf_measurement_model<Measurement, Scalar, NB, NY>
class mekf
{
    static_assert(NB >= 3, "Bias dimension NB must be at least 3: the propagation subtracts the leading three bias elements from the gyro rate");
    static_assert(NY > 0, "Output dimension NY must be positive");

    static constexpr std::size_t NE = 3 + NB;
    static constexpr int ne = static_cast<int>(NE);
    static constexpr int ny = static_cast<int>(NY);
    static constexpr int nb = static_cast<int>(NB);

public:
    using observer_tag = struct mekf_tag;
    using state_vector_t = Vector<Scalar, 4 + NB>;
    using input_vector_t = Vector<Scalar, 3>;
    using output_vector_t = Vector<Scalar, NY>;
    using cov_matrix_t = Matrix<Scalar, NE, NE>;
    using meas_cov_t = Matrix<Scalar, NY, NY>;

    /// @brief Fallible factory. Validates the initial quaternion before the
    /// normalization that seeds the filter state.
    ///
    /// A zero or non-finite `config.q0` norm makes `q0.normalized()` produce
    /// NaN, which would silently poison the whole filter state at construction;
    /// such a config is rejected with `filter_error::degenerate_quaternion`.
    /// Any finite nonzero quaternion is accepted and normalized.
    static auto create(Measurement measurement, mekf_config<Scalar, NB, NY> config) -> ctrlpp::expected<mekf, filter_error>
    {
        const Scalar q0_norm = config.q0.norm();
        if(!(q0_norm > Scalar{0}) || !std::isfinite(q0_norm))
            return ctrlpp::unexpected(filter_error::degenerate_quaternion);
        return mekf{validated_tag{}, std::move(measurement), std::move(config)};
    }

    void predict(const input_vector_t& omega) { predict_impl(omega, dt_); }

    void predict(const input_vector_t& omega, Scalar dt) { predict_impl(omega, dt); }

    /// @brief Update the nominal state and the error-state covariance with a
    /// measurement.
    ///
    /// The step is rejected before any member is assigned when the carried
    /// estimate or the measurement is non-finite, so a rejected step leaves the
    /// attitude, the bias, the covariance and the innovation bitwise unchanged
    /// and the caller may retry with the next sample.
    ///
    /// `predict` is deliberately not fallible: its input is a gyro rate the
    /// caller already owns, and rejecting it would leave the filter with no
    /// propagation for a rotation the body did perform. A prediction that
    /// poisons the nominal state is instead reported by `health()`, which the
    /// next update latches.
    auto update(const output_vector_t& z) -> ctrlpp::expected<void, mekf_update_error>
    {
        if(const auto step = check_step(z); !step)
            return ctrlpp::unexpected(latch_health(step.error()));

        auto z_pred = measurement_(q_, b_);
        innovation_ = (z - z_pred).eval();

        auto H = compute_measurement_jacobian();
        auto S = compute_innovation_covariance(H);
        auto K = compute_kalman_gain(H, S);

        auto delta_xi = apply_multiplicative_correction(K);
        update_covariance(K, H, delta_xi);
        update_state_cache();
        return {};
    }

    auto state() const -> const state_vector_t& { return state_cache_; }
    auto covariance() const -> const cov_matrix_t& { return P_; }
    auto innovation() const -> const output_vector_t& { return innovation_; }
    auto attitude() const -> Eigen::Quaternion<Scalar> { return q_; }
    auto bias() const -> const Vector<Scalar, NB>& { return b_; }

    /// @brief Report whether the carried estimate is still degraded from an
    /// earlier step. Latches; a rejected measurement does not set it.
    auto health() const -> mekf_health { return health_; }

private:
    struct validated_tag
    {
    };

    /// @brief Classify a step's operands without touching a single member.
    ///
    /// The order is the severity order documented on `mekf_update_error`: the
    /// carried estimate first, the supplied measurement last. The cost is one
    /// finiteness scan of each operand -- 4 + NB + NE*NE + NY reads, no branches
    /// on data and no allocation, all dimensions being compile-time constants.
    auto check_step(const output_vector_t& z) const -> ctrlpp::expected<void, mekf_update_error>
    {
        if(!q_.coeffs().allFinite() || !b_.allFinite())
            return ctrlpp::unexpected(mekf_update_error::non_finite_state);
        if(!P_.allFinite())
            return ctrlpp::unexpected(mekf_update_error::non_finite_covariance);
        if(!z.allFinite())
            return ctrlpp::unexpected(mekf_update_error::non_finite_measurement);
        return {};
    }

    /// @brief Latch the persistent status for the faults that describe the
    /// carried estimate, and pass the fault through so the caller still receives
    /// the specific cause on the failure channel.
    auto latch_health(mekf_update_error fault) -> mekf_update_error
    {
        if(fault != mekf_update_error::non_finite_measurement)
            health_ = mekf_health::non_finite_estimate;
        return fault;
    }

    mekf(validated_tag, Measurement measurement, mekf_config<Scalar, NB, NY> config)
        : measurement_{std::move(measurement)}
        , q_{config.q0.normalized()}
        , b_{std::move(config.b0)}
        , P_{std::move(config.P0)}
        , Q_{std::move(config.Q)}
        , R_{std::move(config.R)}
        , eps_{config.numerical_eps}
        , dt_{config.dt}
        , innovation_{output_vector_t::Zero()}
        , state_cache_{}
    {
        update_state_cache();
    }

    /// @brief Propagate nominal quaternion and error-state covariance.
    ///
    /// @cite markley2003 -- Markley, "Attitude Error Representations for Kalman Filtering", 2003
    /// @cite sola2018 -- Sola et al., "A micro Lie theory for state estimation in robotics", 2018
    void predict_impl(const input_vector_t& omega, Scalar dt)
    {
        Vector<Scalar, 3> omega_corr = omega - b_.template head<3>();
        Vector<Scalar, 3> omega_dt = (omega_corr * dt).eval();
        q_ = (q_ * so3::exp(omega_dt)).normalized();

        Eigen::Matrix<Scalar, 3, 3> C = so3::exp(omega_dt).toRotationMatrix();

        // The filter carries a right (body-frame) multiplicative error
        // (q = q_nominal * exp(delta_att), see apply_multiplicative_correction),
        // so the attitude error propagates with the transpose of the incremental
        // rotation: from exp(-omega_dt) * exp(delta) * exp(omega_dt) = exp(C^T delta)
        // the attitude sub-block of F is C^T, not C.
        cov_matrix_t F = cov_matrix_t::Identity();
        F.template block<3, 3>(0, 0) = C.transpose();
        F.template block<3, nb>(0, 3) = -Eigen::Matrix<Scalar, 3, nb>::Identity() * dt;

        P_ = detail::symmetrize((F * P_ * F.transpose() + Q_).eval());
        update_state_cache();
    }

    /// @brief Compute measurement Jacobian H (analytical or numerical).
    ///
    /// @cite markley2003 -- Markley, "Attitude Error Representations for Kalman Filtering", 2003, Eq. 34
    auto compute_measurement_jacobian() const -> Matrix<Scalar, NY, NE>
    {
        if constexpr(differentiable_mekf_measurement<Measurement, Scalar, NB, NY>)
            return measurement_.jacobian(q_, b_);
        else
            return detail::numerical_mekf_jacobian<Scalar, NB, NY>(measurement_, q_, b_, eps_);
    }

    /// @brief Compute innovation covariance: S = H*P*H^T + R.
    auto compute_innovation_covariance(const Matrix<Scalar, NY, NE>& H) const -> meas_cov_t
    {
        return (H * P_ * H.transpose() + R_).eval();
    }

    /// @brief Compute Kalman gain via transpose-solve.
    ///
    /// @cite markley2003 -- Markley, "Attitude Error Representations for Kalman Filtering", 2003
    auto compute_kalman_gain(const Matrix<Scalar, NY, NE>& H, const meas_cov_t& S) const -> Eigen::Matrix<Scalar, ne, ny>
    {
        Eigen::Matrix<Scalar, ny, ne> KT = S.transpose().colPivHouseholderQr().solve(H * P_);
        return KT.transpose().eval();
    }

    /// @brief Apply multiplicative quaternion correction and bias update.
    ///
    /// @cite markley2003 -- Markley, "Attitude Error Representations for Kalman Filtering", 2003, Eq. 46
    auto apply_multiplicative_correction(const Eigen::Matrix<Scalar, ne, ny>& K) -> Vector<Scalar, NE>
    {
        Vector<Scalar, NE> delta_xi = K * innovation_;
        Vector<Scalar, 3> delta_att = delta_xi.template head<3>().eval();
        q_ = (q_ * so3::exp(delta_att)).normalized();
        b_ += delta_xi.template tail<NB>();
        return delta_xi;
    }

    /// @brief Joseph-form covariance update with mandatory frame-change reset.
    ///
    /// @cite markley2003 -- Markley, "Attitude Error Representations for Kalman Filtering", 2003, Eq. 68
    void update_covariance(const Eigen::Matrix<Scalar, ne, ny>& K, const Matrix<Scalar, NY, NE>& H, const Vector<Scalar, NE>& delta_xi)
    {
        cov_matrix_t IKH = cov_matrix_t::Identity() - K * H;
        P_ = (IKH * P_ * IKH.transpose() + K * R_ * K.transpose()).eval();

        // Mandatory covariance reset via frame-change Jacobian G
        Vector<Scalar, 3> delta_att = delta_xi.template head<3>();
        cov_matrix_t G = cov_matrix_t::Identity();
        G.template block<3, 3>(0, 0) -= Scalar{0.5} * so3::skew(delta_att);
        P_ = detail::symmetrize((G * P_ * G.transpose()).eval());
    }

    void update_state_cache()
    {
        state_cache_(0) = q_.w();
        state_cache_(1) = q_.x();
        state_cache_(2) = q_.y();
        state_cache_(3) = q_.z();
        for(std::size_t i = 0; i < NB; ++i)
        {
            state_cache_(static_cast<Eigen::Index>(4 + i)) = b_(static_cast<Eigen::Index>(i));
        }
    }

    Measurement measurement_;
    Eigen::Quaternion<Scalar> q_;
    Vector<Scalar, NB> b_;
    cov_matrix_t P_;
    cov_matrix_t Q_;
    meas_cov_t R_;
    Scalar eps_;
    Scalar dt_;
    output_vector_t innovation_;
    state_vector_t state_cache_;
    mekf_health health_{mekf_health::ok};
};

namespace detail
{

struct mekf_sa_measurement
{
    auto operator()(const Eigen::Quaternion<double>&, const Vector<double, 3>&) const -> Vector<double, 3> { return Vector<double, 3>::Zero(); }
};

}

static_assert(ObserverPolicy<mekf<double, 3, 3, detail::mekf_sa_measurement>>);
static_assert(CovarianceObserver<mekf<double, 3, 3, detail::mekf_sa_measurement>>);

}

#endif
