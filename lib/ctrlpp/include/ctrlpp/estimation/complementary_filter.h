#ifndef HPP_GUARD_CTRLPP_ESTIMATION_COMPLEMENTARY_FILTER_H
#define HPP_GUARD_CTRLPP_ESTIMATION_COMPLEMENTARY_FILTER_H

/// @brief Mahony nonlinear complementary filter for attitude estimation on SO(3).
///
/// Provides computationally lightweight attitude estimation from IMU (gyro + accel)
/// or MARG (gyro + accel + mag) sensor data. Satisfies ObserverPolicy.
///
/// @cite mahony2008 -- Mahony et al., "Nonlinear Complementary Filters on the Special Orthogonal Group", 2008

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/lie/so3.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/estimation_types.h"

#include <Eigen/Geometry>

#include <cmath>
#include <utility>

namespace ctrlpp
{

/// @brief Structured failure modes of a `complementary_filter` update.
///
/// Each enumerator is an exact domain condition, not a tuning preference: a
/// non-finite operand makes every downstream product non-finite, so the step
/// cannot produce an attitude at all. The filter carries no covariance, so the
/// carried estimate is the attitude quaternion together with the gyro bias.
///
///  * non_finite_state       : the carried attitude quaternion or the carried
///                             gyro bias is already non-finite when the step
///                             begins. It is reported ahead of the sensor
///                             vectors because the fault is upstream of them,
///                             and because both correction terms are computed
///                             from the carried attitude's rotation matrix.
///  * non_finite_measurement : one of the supplied sensor vectors -- rate,
///                             acceleration, or magnetic field where present --
///                             has a non-finite component. Every one of them
///                             reaches the quaternion integration, whose output
///                             is the filter's carried memory, so a single such
///                             sample destroys the attitude permanently.
///  * non_finite_timestep    : the supplied integration step is not finite. A
///                             distinguishable cause from a non-finite sensor
///                             reading because it names a broken clock rather
///                             than a broken sensor, and the caller fixes a
///                             different input. It poisons the integration just
///                             as surely: the tangent vector is the step times
///                             the corrected rate.
enum class cf_update_error
{
    non_finite_state,
    non_finite_measurement,
    non_finite_timestep,
    non_finite_result,
};

/// @brief Persistent state-health status of a `complementary_filter`.
///
/// A per-call result cannot answer whether the carried attitude is still
/// degraded from a step several samples ago, because that question outlives the
/// call. The status latches.
///
///  * ok                  : every step so far began from a finite attitude and
///                          bias.
///  * non_finite_estimate : a step found the carried attitude or bias already
///                          non-finite. A rejected update does NOT set it: the
///                          rejection mutates nothing, so it leaves the filter
///                          healthy.
enum class cf_health
{
    ok,
    non_finite_estimate,
};

template <ctrlpp_floating_scalar Scalar>
struct cf_config
{
    Scalar k_p{Scalar{2}};
    Scalar k_i{static_cast<Scalar>(0.005)};
    Scalar dt{static_cast<Scalar>(0.01)};
    Eigen::Quaternion<Scalar> q0{Eigen::Quaternion<Scalar>::Identity()};
};

template <ctrlpp_floating_scalar Scalar>
class complementary_filter
{
public:
    using observer_tag = struct complementary_filter_tag;
    using state_vector_t = Vector<Scalar, 7>;
    using input_vector_t = Vector<Scalar, 3>;
    using output_vector_t = Vector<Scalar, 3>;

    /// @brief Fallible factory. Validates and normalizes the initial
    /// quaternion before it seeds the filter state.
    ///
    /// A zero or non-finite `config.q0` norm cannot be normalized without
    /// producing NaN, so such a config is rejected with
    /// `filter_error::degenerate_quaternion`. Any finite nonzero quaternion is
    /// accepted and normalized: the correction terms treat the stored
    /// quaternion as a unit rotation via `toRotationMatrix()`, so a non-unit
    /// q0 is brought onto the unit sphere at construction (previously it was
    /// stored raw and a non-unit q0 skewed the first gravity and magnetic
    /// references until the first gyro integration renormalized it).
    static auto create(cf_config<Scalar> config) -> ctrlpp::expected<complementary_filter, filter_error>
    {
        const Scalar q0_norm = config.q0.norm();
        if(!(q0_norm > Scalar{0}) || !std::isfinite(q0_norm))
            return ctrlpp::unexpected(filter_error::degenerate_quaternion);
        return complementary_filter{validated_tag{}, std::move(config)};
    }

    // Natural IMU update (6-DOF): gyro + accelerometer.
    ///
    /// The step is rejected before any member is assigned when the carried
    /// attitude, either sensor vector, or the integration step is non-finite, so
    /// a rejected step leaves the attitude and the bias bitwise unchanged and
    /// the caller may retry with the next sample.
    ///
    /// A directionless acceleration is NOT a rejection and NOT a skipped step.
    /// The step succeeded: a zero acceleration vector carries no gravity
    /// direction, so there is no correction to apply, and the rate is
    /// integrated with a zero correction exactly as the algorithm prescribes.
    /// The body rotated whether or not the accelerometer could say which way is
    /// down, so discarding the integration would drop a rotation that happened
    /// and report the step a success. Reporting the missing correction through
    /// the failure channel would be the opposite error: it would tell the
    /// caller their step failed when it did exactly what the algorithm
    /// prescribes, and a caller who learns the failure channel carries
    /// non-failures will eventually ignore a real one.
    ///
    /// @cite mahony2008 -- Mahony et al., 2008, Sec. III (IMU complementary filter)
    auto update(const Vector<Scalar, 3>& gyro, const Vector<Scalar, 3>& accel, Scalar dt) -> ctrlpp::expected<void, cf_update_error>
    {
        if(const auto step = check_step(gyro, accel, dt); !step)
            return ctrlpp::unexpected(latch_health(step.error()));

        auto const previous_q = q_;
        auto const previous_bias = bias_;
        auto const previous_state = state_cache_;
        Vector<Scalar, 3> direction;
        Vector<Scalar, 3> correction = Vector<Scalar, 3>::Zero();
        if(unit_direction(accel, direction))
            correction = compute_gravity_correction(direction);

        integrate_gyro(gyro, correction, dt);
        if(!q_.coeffs().allFinite() || !bias_.allFinite()
            || !state_cache_.allFinite())
        {
            q_ = previous_q;
            bias_ = previous_bias;
            state_cache_ = previous_state;
            return ctrlpp::unexpected(cf_update_error::non_finite_result);
        }
        return {};
    }

    // Natural MARG update (9-DOF): gyro + accelerometer + magnetometer.
    ///
    /// Rejects on the same terms as the IMU overload, with the magnetic vector
    /// added to the sensor operands. The two corrections are independent sums in
    /// Mahony's law, so each is applied when its own vector has a direction and
    /// omitted when it does not, and the rate is integrated either way. That is
    /// the same shape as the six-axis overload with one correction instead of
    /// two.
    ///
    /// @cite mahony2008 -- Mahony et al., 2008, Sec. IV (MARG complementary filter)
    auto update(const Vector<Scalar, 3>& gyro, const Vector<Scalar, 3>& accel, const Vector<Scalar, 3>& mag, Scalar dt) -> ctrlpp::expected<void, cf_update_error>
    {
        if(const auto step = check_step(gyro, accel, dt, &mag); !step)
            return ctrlpp::unexpected(latch_health(step.error()));

        auto const previous_q = q_;
        auto const previous_bias = bias_;
        auto const previous_state = state_cache_;
        Vector<Scalar, 3> direction;
        Vector<Scalar, 3> correction = Vector<Scalar, 3>::Zero();
        if(unit_direction(accel, direction))
            correction += compute_gravity_correction(direction);
        if(unit_direction(mag, direction))
            correction += compute_magnetic_correction(direction);

        integrate_gyro(gyro, correction, dt);
        if(!q_.coeffs().allFinite() || !bias_.allFinite()
            || !state_cache_.allFinite())
        {
            q_ = previous_q;
            bias_ = previous_bias;
            state_cache_ = previous_state;
            return ctrlpp::unexpected(cf_update_error::non_finite_result);
        }
        return {};
    }

    // ObserverPolicy wrappers (use config dt)
    void predict(const input_vector_t& u) { gyro_buf_ = u; }

    /// @brief Observer-concept form: correct with the buffered rate from the
    /// last `predict` and the configured step. Forwards the IMU overload's
    /// result rather than swallowing it, so the caller sees the same cause.
    auto update(const output_vector_t& z) -> ctrlpp::expected<void, cf_update_error> { return update(gyro_buf_, z, dt_); }

    auto state() const -> const state_vector_t& { return state_cache_; }
    auto attitude() const -> Eigen::Quaternion<Scalar> { return q_; }
    auto bias() const -> const Vector<Scalar, 3>& { return bias_; }

    /// @brief Report whether the carried attitude is still degraded from an
    /// earlier step. Latches; a rejected update does not set it.
    auto health() const -> cf_health { return health_; }

private:
    struct validated_tag
    {
    };

    /// @brief Classify a step's operands without touching a single member.
    ///
    /// The order is the severity order documented on `cf_update_error`: the
    /// carried attitude and bias first, then the sensor vectors, then the
    /// integration step. The magnetic vector is checked only when the caller
    /// supplied one, which the trailing parameter's pointer encodes. The cost is
    /// one finiteness scan of each operand -- at most 7 + 9 + 1 reads, no
    /// branches on data and no allocation.
    auto check_step(const Vector<Scalar, 3>& gyro, const Vector<Scalar, 3>& accel, Scalar dt, const Vector<Scalar, 3>* mag = nullptr) const -> ctrlpp::expected<void, cf_update_error>
    {
        if(!q_.coeffs().allFinite() || !bias_.allFinite())
            return ctrlpp::unexpected(cf_update_error::non_finite_state);
        if(!gyro.allFinite() || !accel.allFinite() || (mag != nullptr && !mag->allFinite()))
            return ctrlpp::unexpected(cf_update_error::non_finite_measurement);
        if(!std::isfinite(dt))
            return ctrlpp::unexpected(cf_update_error::non_finite_timestep);
        return {};
    }

    /// @brief Write the unit direction of a sensor reading and report whether it
    /// has one.
    ///
    /// The coefficients are divided by their largest magnitude before the norm
    /// is formed, so the norm is taken of a vector whose largest coefficient is
    /// exactly one and whose squared norm therefore lies in [1, 3]. Neither end
    /// of the exponent range is reachable from there, so every finite nonzero
    /// reading has a direction -- including one whose squared norm would
    /// underflow to zero, or overflow to infinity, if the norm were formed
    /// directly from the reading. The only reading with no direction is the
    /// exactly zero vector, which carries none at all rather than a small one.
    ///
    /// There is deliberately NO magnitude threshold here. A sensor reading is
    /// expressed in the caller's units, so an absolute floor gives the same
    /// physical acceleration different treatment depending on whether it is
    /// reported in g or in millimeters per second squared -- and both
    /// correction terms use only the DIRECTION, which is scale-free. The
    /// previous floor of 1e-10 discarded a reading of 1e-15 g, which is an
    /// ordinary reading in units chosen that way.
    ///
    /// The returned bool is a predicate on the argument, not a failure channel:
    /// the only reason a vector has no direction is that every coefficient is
    /// zero, so there is no cause to carry and nothing is lost by a bool.
    ///
    /// This mirrors `so3::normalize`, which repairs the identical defect on the
    /// attitude quaternion and for the identical reason.
    static auto unit_direction(const Vector<Scalar, 3>& v, Vector<Scalar, 3>& direction) -> bool
    {
        const Scalar scale = v.cwiseAbs().maxCoeff();
        if(!(scale > Scalar{0}))
            return false;
        const Vector<Scalar, 3> scaled = (v / scale).eval();
        direction = (scaled / scaled.norm()).eval();
        return true;
    }

    /// @brief Latch the persistent status for the fault that describes the
    /// carried attitude, and pass the fault through so the caller still receives
    /// the specific cause on the failure channel.
    auto latch_health(cf_update_error fault) -> cf_update_error
    {
        if(fault == cf_update_error::non_finite_state)
            health_ = cf_health::non_finite_estimate;
        return fault;
    }

    complementary_filter(validated_tag, cf_config<Scalar> config) : q_{config.q0.normalized()}, bias_{Vector<Scalar, 3>::Zero()}, k_p_{config.k_p}, k_i_{config.k_i}, dt_{config.dt}, gyro_buf_{Vector<Scalar, 3>::Zero()}
    {
        update_state_cache();
    }

    /// @cite mahony2008 -- Mahony et al., 2008, Eq. 12 (gravity error via cross product)
    auto compute_gravity_correction(const Vector<Scalar, 3>& acc_n) const -> Vector<Scalar, 3>
    {
        Vector<Scalar, 3> g_body = q_.toRotationMatrix().transpose().col(2);
        return acc_n.cross(g_body);
    }

    /// @cite mahony2008 -- Mahony et al., 2008, Eq. 48 (magnetic field error)
    auto compute_magnetic_correction(const Vector<Scalar, 3>& mag_n) const -> Vector<Scalar, 3>
    {
        auto R = q_.toRotationMatrix();
        Vector<Scalar, 3> h = R * mag_n;
        Vector<Scalar, 3> b{std::sqrt(h(0) * h(0) + h(1) * h(1)), Scalar{0}, h(2)};
        Vector<Scalar, 3> b_body = R.transpose() * b;
        return mag_n.cross(b_body);
    }

    /// @cite mahony2008 -- Mahony et al., 2008, Eq. 6 (quaternion integration with PI correction)
    void integrate_gyro(const Vector<Scalar, 3>& gyro, const Vector<Scalar, 3>& e, Scalar dt)
    {
        bias_ -= k_i_ * e * dt;
        Vector<Scalar, 3> omega_c = gyro - bias_ + k_p_ * e;
        Vector<Scalar, 3> phi = dt * omega_c;
        q_ = (q_ * so3::exp(phi)).normalized();
        update_state_cache();
    }

    void update_state_cache() { state_cache_ << q_.w(), q_.x(), q_.y(), q_.z(), bias_(0), bias_(1), bias_(2); }

    Eigen::Quaternion<Scalar> q_;
    Vector<Scalar, 3> bias_;
    Scalar k_p_;
    Scalar k_i_;
    Scalar dt_;
    Vector<Scalar, 3> gyro_buf_;
    state_vector_t state_cache_;
    cf_health health_{cf_health::ok};
};

static_assert(ObserverPolicy<complementary_filter<double>>);

}

#endif
