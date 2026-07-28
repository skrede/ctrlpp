// kalman_tracking.cpp -- Kalman filter tracking matching kalman_tracking.m
// Usage: ./kalman_tracking > kalman_tracking_cpp.csv

#include "ctrlpp/estimation/kalman.h"
#include "ctrlpp/model/discretise.h"
#include "ctrlpp/model/propagate.h"
#include "ctrlpp/model/state_space.h"
#include "ctrlpp/control/lqr.h"

#include <cstdio>
#include <iostream>

int main()
{
    using Scalar = double;
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;
    constexpr std::size_t NY = 1;

    ctrlpp::continuous_state_space<Scalar, NX, NU, NY> sys_c{};
    sys_c.A << 0.0, 1.0, -1.0, -0.5;
    sys_c.B << 0.0, 1.0;
    sys_c.C << 1.0, 0.0;
    sys_c.D << 0.0;

    constexpr Scalar dt = 0.05;
    constexpr Scalar duration = 10.0;

    auto sys_d = ctrlpp::discretise(ctrlpp::zoh{}, sys_c, dt);

    Eigen::Matrix<Scalar, 2, 2> Q_proc = Eigen::Matrix<Scalar, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<Scalar, 1, 1> R_meas;
    R_meas << 0.1;
    Eigen::Matrix<Scalar, 2, 2> P0 = Eigen::Matrix<Scalar, 2, 2>::Identity();
    Eigen::Matrix<Scalar, 2, 1> x0_est = Eigen::Matrix<Scalar, 2, 1>::Zero();

    auto kf_result = ctrlpp::kalman_filter<Scalar, NX, NU, NY>::create(sys_d, {.Q = Q_proc, .R = R_meas, .x0 = x0_est, .P0 = P0});
    if(!kf_result.has_value())
    {
        std::cerr << "invalid Kalman filter configuration\n";
        return 1;
    }
    auto& kf = *kf_result;

    Eigen::Matrix<Scalar, 2, 2> Q_lqr;
    Q_lqr << 10.0, 0.0, 0.0, 1.0;
    Eigen::Matrix<Scalar, 1, 1> R_lqr;
    R_lqr << 1.0;
    auto K_result = ctrlpp::lqr_gain<Scalar, NX, NU>(sys_d.A, sys_d.B, Q_lqr, R_lqr);
    if(!K_result.has_value())
    {
        std::cerr << "LQR gain synthesis declined the validation plant\n";
        return 1;
    }
    ctrlpp::lqr<Scalar, NX, NU> controller(*K_result);

    Eigen::Matrix<Scalar, 2, 1> x_true;
    x_true << 1.0, 0.0;

    std::printf("time,x_true_0,x_true_1,x_est_0,x_est_1,P_00,P_11\n");

    for(Scalar t = 0.0; t < duration - dt / 2; t += dt)
    {
        auto x_est = kf.state();
        auto u = controller.compute(x_est);
        Eigen::Matrix<Scalar, 1, 1> z = sys_d.C * x_true;

        auto P = kf.covariance();
        std::printf("%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n",
                    t, x_true(0), x_true(1), x_est(0), x_est(1), P(0, 0), P(1, 1));

        kf.predict(u);
        x_true = ctrlpp::propagate(sys_d, x_true, u);
        Eigen::Matrix<Scalar, 1, 1> z_new = sys_d.C * x_true;
        if(!kf.update(z_new).has_value())
        {
            std::cerr << "the filter refused a sample\n";
            return 1;
        }
    }
}
