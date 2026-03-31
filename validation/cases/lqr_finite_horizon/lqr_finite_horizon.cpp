// lqr_finite_horizon.cpp -- Finite-horizon LQR matching lqr_finite_horizon.m
// Usage: ./lqr_finite_horizon > lqr_finite_horizon_cpp.csv

#include "ctrlpp/model/discretise.h"
#include "ctrlpp/model/propagate.h"
#include "ctrlpp/model/state_space.h"
#include "ctrlpp/control/lqr.h"

#include <cstdio>

int main()
{
    using Scalar = double;
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;
    constexpr std::size_t NY = 2;

    ctrlpp::continuous_state_space<Scalar, NX, NU, NY> sys_c{};
    sys_c.A << 0.0, 1.0, -1.0, -0.5;
    sys_c.B << 0.0, 1.0;
    sys_c.C = Eigen::Matrix<Scalar, 2, 2>::Identity();
    sys_c.D = Eigen::Matrix<Scalar, 2, 1>::Zero();

    constexpr Scalar dt = 0.05;

    auto sys_d = ctrlpp::discretise(ctrlpp::zoh{}, sys_c, dt);

    Eigen::Matrix<Scalar, 2, 2> Q;
    Q << 10.0, 0.0, 0.0, 1.0;
    Eigen::Matrix<Scalar, 1, 1> R;
    R << 1.0;
    Eigen::Matrix<Scalar, 2, 2> Qf;
    Qf << 20.0, 0.0, 0.0, 2.0;
    constexpr std::size_t horizon = 50;

    auto gains = ctrlpp::lqr_finite<Scalar, NX, NU>(sys_d.A, sys_d.B, Q, R, Qf, horizon);

    Eigen::Matrix<Scalar, 2, 1> x;
    x << 1.0, 0.0;

    std::printf("step,x0,x1,u,K0,K1\n");

    for(std::size_t k = 0; k < horizon; ++k)
    {
        auto u = (-gains[k] * x).eval();
        std::printf("%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n",
                    static_cast<Scalar>(k + 1), x(0), x(1), u(0),
                    gains[k](0, 0), gains[k](0, 1));
        Eigen::Matrix<Scalar, 1, 1> u_vec;
        u_vec << u(0);
        x = ctrlpp::propagate(sys_d, x, u_vec);
    }
}
