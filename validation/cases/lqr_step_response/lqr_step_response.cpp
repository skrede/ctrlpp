// lqr_step_response.cpp -- LQR closed-loop step response matching lqr_step_response.m
// Usage: ./lqr_step_response > lqr_step_response_cpp.csv

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
    constexpr Scalar duration = 10.0;

    auto sys_d = ctrlpp::discretise(ctrlpp::zoh{}, sys_c, dt);

    Eigen::Matrix<Scalar, 2, 2> Q;
    Q << 10.0, 0.0, 0.0, 1.0;
    Eigen::Matrix<Scalar, 1, 1> R;
    R << 1.0;

    ctrlpp::lqr<Scalar, NX, NU> controller(*ctrlpp::lqr_gain<Scalar, NX, NU>(sys_d.A, sys_d.B, Q, R));

    Eigen::Matrix<Scalar, 2, 1> x;
    x << 1.0, 0.0;

    std::printf("time,x0,x1,u\n");

    for(Scalar t = 0.0; t < duration - dt / 2; t += dt)
    {
        auto u = controller.compute(x);
        std::printf("%.15e,%.15e,%.15e,%.15e\n", t, x(0), x(1), u(0));
        x = ctrlpp::propagate(sys_d, x, u);
    }
}
