// Embedded-clean compile witness: instantiates the core ctrlpp surface with
// Scalar=float and exercises one representative create per converted module
// family. Every fallible result is handled through has_value() plus operator*,
// so the unit stays free of the exception machinery and compiles under
// -fno-exceptions with CTRLPP_NO_EXCEPTIONS defined. It is registered in the
// host build so the embeddable surface cannot silently rot.

#include "ctrlpp/control.h"
#include "ctrlpp/estimation.h"
#include "ctrlpp/dsp.h"
#include "ctrlpp/trajectory.h"
#include "ctrlpp/model/state_space.h"

#include <Eigen/Dense>

#include <random>
#include <cstddef>

namespace
{

using scalar = float;

// -- Minimal Euclidean models satisfying the estimator concepts --

struct linear_dynamics
{
    auto operator()(const ctrlpp::Vector<scalar, 2>& x, const ctrlpp::Vector<scalar, 1>& u) const
        -> ctrlpp::Vector<scalar, 2>
    {
        ctrlpp::Vector<scalar, 2> next;
        next[0] = x[0] + u[0];
        next[1] = x[1];
        return next;
    }

    auto jacobian_x(const ctrlpp::Vector<scalar, 2>&, const ctrlpp::Vector<scalar, 1>&) const
        -> ctrlpp::Matrix<scalar, 2, 2>
    {
        return ctrlpp::Matrix<scalar, 2, 2>::Identity();
    }

    auto jacobian_u(const ctrlpp::Vector<scalar, 2>&, const ctrlpp::Vector<scalar, 1>&) const
        -> ctrlpp::Matrix<scalar, 2, 1>
    {
        return ctrlpp::Matrix<scalar, 2, 1>::Zero();
    }
};

struct linear_measurement
{
    auto operator()(const ctrlpp::Vector<scalar, 2>& x) const -> ctrlpp::Vector<scalar, 1>
    {
        return ctrlpp::Vector<scalar, 1>{x[0]};
    }

    auto jacobian(const ctrlpp::Vector<scalar, 2>&) const -> ctrlpp::Matrix<scalar, 1, 2>
    {
        ctrlpp::Matrix<scalar, 1, 2> h;
        h << scalar{1}, scalar{0};
        return h;
    }
};

// -- Minimal manifold models for the SO(3) estimators --

struct attitude_dynamics
{
    auto operator()(const Eigen::Quaternion<scalar>& q, const ctrlpp::Vector<scalar, 3>&) const
        -> Eigen::Quaternion<scalar>
    {
        return q;
    }
};

struct attitude_measurement
{
    auto operator()(const Eigen::Quaternion<scalar>&) const -> ctrlpp::Vector<scalar, 3>
    {
        return ctrlpp::Vector<scalar, 3>::Zero();
    }
};

struct mekf_gravity_measurement
{
    auto operator()(const Eigen::Quaternion<scalar>&, const ctrlpp::Vector<scalar, 3>&) const
        -> ctrlpp::Vector<scalar, 3>
    {
        return ctrlpp::Vector<scalar, 3>::Zero();
    }
};

// -- Fold any expected/optional into an integer so nothing is optimized away --

template <typename Expected>
auto fold(const Expected& result) -> int
{
    return result.has_value() ? 1 : 0;
}

}

int main()
{
    int witness = 0;

    // -- control: pid, lqr gain, dare, care --

    ctrlpp::pid_config<scalar, 1> pid_cfg{};
    pid_cfg.kp = ctrlpp::Vector<scalar, 1>::Constant(scalar{1});
    ctrlpp::pid<scalar, 1> controller(pid_cfg);
    const auto u = controller.compute(
        ctrlpp::Vector<scalar, 1>::Zero(), ctrlpp::Vector<scalar, 1>::Zero(), scalar{0.01});
    witness += (u.has_value() && (*u)[0] == (*u)[0]) ? 1 : 0;

    const ctrlpp::Matrix<scalar, 2, 2> a = ctrlpp::Matrix<scalar, 2, 2>::Identity();
    ctrlpp::Matrix<scalar, 2, 1> b;
    b << scalar{0}, scalar{1};
    const ctrlpp::Matrix<scalar, 2, 2> q = ctrlpp::Matrix<scalar, 2, 2>::Identity();
    ctrlpp::Matrix<scalar, 1, 1> r;
    r << scalar{1};

    const auto k = ctrlpp::lqr_gain<scalar, 2, 1>(a, b, q, r);
    witness += k.has_value() ? 1 : 0;

    witness += fold(ctrlpp::dare<scalar, 2, 1>(a, b, q, r));
    witness += fold(ctrlpp::care<scalar, 2, 1>(a, b, q, r));

    // -- estimation: plain-ctor filters (particle) --

    auto pf = ctrlpp::make_particle_filter<8>(
        linear_dynamics{}, linear_measurement{}, ctrlpp::pf_config<scalar, 2, 1, 1>{},
        std::mt19937_64{42U});
    witness += static_cast<int>(sizeof(pf) > 0);

    // -- estimation: create factories (kalman, ekf, ukf, mekf, manifold_ukf, complementary) --

    ctrlpp::discrete_state_space<scalar, 2, 1, 1> system{};
    witness += fold(
        ctrlpp::kalman_filter<scalar, 2, 1, 1>::create(system, ctrlpp::kalman_config<scalar, 2, 1, 1>{}));

    witness += fold(
        ctrlpp::ekf<scalar, 2, 1, 1, linear_dynamics, linear_measurement>::create(
            linear_dynamics{}, linear_measurement{}, ctrlpp::ekf_config<scalar, 2, 1, 1>{}));

    witness += fold(
        ctrlpp::ukf<scalar, 2, 1, 1, linear_dynamics, linear_measurement>::create(
            linear_dynamics{}, linear_measurement{}, ctrlpp::ukf_config<scalar, 2, 1, 1>{}));

    witness += fold(
        ctrlpp::mekf<scalar, 3, 3, mekf_gravity_measurement>::create(
            mekf_gravity_measurement{}, ctrlpp::mekf_config<scalar, 3, 3>{}));

    witness += fold(
        ctrlpp::manifold_ukf<scalar, 3, attitude_dynamics, attitude_measurement>::create(
            attitude_dynamics{}, attitude_measurement{}, ctrlpp::manifold_ukf_config<scalar, 3>{}));

    witness += fold(ctrlpp::complementary_filter<scalar>::create(ctrlpp::cf_config<scalar>{}));

    // -- dsp: biquad expected factory + fir plain ctor --

    witness += fold(ctrlpp::biquad<scalar>::low_pass(scalar{100}, scalar{1000}));

    const ctrlpp::fir<scalar, 3> filter({scalar{0.25}, scalar{0.5}, scalar{0.25}});
    witness += static_cast<int>(sizeof(filter) > 0);

    // -- trajectory: cubic/trapezoidal/double_s expected factories, spline + planner create --

    const ctrlpp::Vector<scalar, 1> q0 = ctrlpp::Vector<scalar, 1>::Zero();
    const ctrlpp::Vector<scalar, 1> q1 = ctrlpp::Vector<scalar, 1>::Constant(scalar{1});
    const ctrlpp::Vector<scalar, 1> v_zero = ctrlpp::Vector<scalar, 1>::Zero();
    witness += fold(ctrlpp::make_cubic_trajectory<scalar, 1>(q0, q1, v_zero, v_zero, scalar{1}));

    witness += fold(ctrlpp::trapezoidal_trajectory<scalar>::create(
        {.q0 = scalar{0}, .q1 = scalar{1}, .v_max = scalar{1}, .a_max = scalar{1}}));

    witness += fold(ctrlpp::double_s_trajectory<scalar>::create(
        {.q0 = scalar{0}, .q1 = scalar{1}, .v_max = scalar{1}, .a_max = scalar{1}, .j_max = scalar{1}}));

    witness += fold(ctrlpp::cubic_spline<scalar>::create(
        {.times = {scalar{0}, scalar{1}, scalar{2}}, .positions = {scalar{0}, scalar{1}, scalar{0}}}));

    witness += fold(ctrlpp::online_planner_3rd<scalar>::create(
        {.v_max = scalar{1}, .a_max = scalar{1}, .j_max = scalar{1}}));

    return witness > 0 ? 0 : 1;
}
