// lqi_tracking.cpp -- LQI tracking matching lqi_tracking.m
// Usage: ./lqi_tracking > lqi_tracking_cpp.csv

#include "ctrlpp/model/discretize.h"
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
    constexpr Scalar duration = 15.0;

    auto sys_d = ctrlpp::discretize(ctrlpp::zoh{}, sys_c, dt);

    Eigen::Matrix<Scalar, 3, 3> Q_aug;
    Q_aug << 10.0, 0.0, 0.0,
             0.0, 1.0, 0.0,
             0.0, 0.0, 50.0;
    Eigen::Matrix<Scalar, 1, 1> R;
    R << 0.1;

    auto lqi_result = ctrlpp::lqi_gain<Scalar, NX, NU, NY>(sys_d.A, sys_d.B, sys_d.C, Q_aug, R);
    if(!lqi_result.has_value())
    {
        std::cerr << "LQI gain synthesis declined the validation plant\n";
        return 1;
    }
    auto lqi_res = *lqi_result;

    Eigen::Matrix<Scalar, 2, 1> x = Eigen::Matrix<Scalar, 2, 1>::Zero();
    Scalar xi = 0.0;
    constexpr Scalar ref = 1.0;

    std::printf("time,position,velocity,integral,control\n");

    for(Scalar t = 0.0; t < duration - dt / 2; t += dt)
    {
        Scalar y = (sys_d.C * x)(0);
        Eigen::Matrix<Scalar, 1, 1> u_vec = -lqi_res.Kx * x;
        Eigen::Matrix<Scalar, 1, 1> xi_vec;
        xi_vec << xi;
        u_vec -= lqi_res.Ki * xi_vec;
        Scalar u = u_vec(0);

        Eigen::Matrix<Scalar, 1, 1> u_input;
        u_input << u;

        std::printf("%.15e,%.15e,%.15e,%.15e,%.15e\n", t, x(0), x(1), xi, u);

        x = ctrlpp::propagate(sys_d, x, u_input);
        xi += ref - y;
    }
}
