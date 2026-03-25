#ifndef HPP_GUARD_CTRLPP_ESTIMATION_COMPLEMENTARY_FILTER_H
#define HPP_GUARD_CTRLPP_ESTIMATION_COMPLEMENTARY_FILTER_H

/// @brief Mahony nonlinear complementary filter for attitude estimation on SO(3).
///
/// Provides computationally lightweight attitude estimation from IMU (gyro + accel)
/// or MARG (gyro + accel + mag) sensor data. Satisfies ObserverPolicy.
///
/// @cite mahony2008 -- Mahony et al., "Nonlinear Complementary Filters on the Special Orthogonal Group", 2008

#include "ctrlpp/lie/so3.h"
#include "ctrlpp/types.h"
#include "ctrlpp/estimation/observer_policy.h"

#include <Eigen/Geometry>

#include <cmath>

namespace ctrlpp
{

template <typename Scalar>
struct cf_config
{
    Scalar k_p{Scalar{2}};
    Scalar k_i{Scalar{0.005}};
    Scalar dt{Scalar{0.01}};
    Eigen::Quaternion<Scalar> q0{Eigen::Quaternion<Scalar>::Identity()};
};

template <typename Scalar>
class complementary_filter
{
public:
    using observer_tag = struct complementary_filter_tag;
    using state_vector_t = Vector<Scalar, 7>;
    using input_vector_t = Vector<Scalar, 3>;
    using output_vector_t = Vector<Scalar, 3>;

    explicit complementary_filter(cf_config<Scalar> config) : q_{config.q0}, bias_{Vector<Scalar, 3>::Zero()}, k_p_{config.k_p}, k_i_{config.k_i}, dt_{config.dt}, gyro_buf_{Vector<Scalar, 3>::Zero()}
    {
        update_state_cache();
    }

    // Natural IMU update (6-DOF): gyro + accelerometer.
    /// @cite mahony2008 -- Mahony et al., 2008, Sec. III (IMU complementary filter)
    void update(const Vector<Scalar, 3>& gyro, const Vector<Scalar, 3>& accel, Scalar dt)
    {
        Scalar norm = accel.norm();
        if(norm < Scalar{1e-10})
            return;

        auto e = compute_gravity_correction(accel / norm);
        integrate_gyro(gyro, e, dt);
    }

    // Natural MARG update (9-DOF): gyro + accelerometer + magnetometer.
    /// @cite mahony2008 -- Mahony et al., 2008, Sec. IV (MARG complementary filter)
    void update(const Vector<Scalar, 3>& gyro, const Vector<Scalar, 3>& accel, const Vector<Scalar, 3>& mag, Scalar dt)
    {
        Scalar norm = accel.norm();
        if(norm < Scalar{1e-10})
            return;

        auto e_acc = compute_gravity_correction(accel / norm);

        Scalar mag_norm = mag.norm();
        if(mag_norm < Scalar{1e-10})
        {
            integrate_gyro(gyro, e_acc, dt);
            return;
        }

        auto e_mag = compute_magnetic_correction(mag / mag_norm);
        integrate_gyro(gyro, e_acc + e_mag, dt);
    }

    // ObserverPolicy wrappers (use config dt)
    void predict(const input_vector_t& u) { gyro_buf_ = u; }
    void update(const output_vector_t& z) { update(gyro_buf_, z, dt_); }

    [[nodiscard]] auto state() const -> const state_vector_t& { return state_cache_; }
    [[nodiscard]] auto attitude() const -> Eigen::Quaternion<Scalar> { return q_; }
    [[nodiscard]] auto bias() const -> const Vector<Scalar, 3>& { return bias_; }

private:
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
        bias_ += k_i_ * e * dt;
        Vector<Scalar, 3> omega_c = gyro - bias_ + k_p_ * e;
        Vector<Scalar, 3> phi = (Scalar{0.5} * dt) * omega_c;
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
};

template <typename Scalar>
complementary_filter(cf_config<Scalar>) -> complementary_filter<Scalar>;

static_assert(ObserverPolicy<complementary_filter<double>>);

} // namespace ctrlpp

#endif
