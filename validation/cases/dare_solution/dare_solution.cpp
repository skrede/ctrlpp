// dare_solution.cpp -- DARE solution matching dare_solution.m
// Usage: ./dare_solution > dare_solution_cpp.csv

#include "ctrlpp/model/discretise.h"
#include "ctrlpp/model/state_space.h"
#include "ctrlpp/control/dare.h"
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

    auto sys_d = ctrlpp::discretise(ctrlpp::zoh{}, sys_c, 0.05);

    Eigen::Matrix<Scalar, 2, 2> Q;
    Q << 10.0, 0.0, 0.0, 1.0;
    Eigen::Matrix<Scalar, 1, 1> R;
    R << 1.0;

    auto P = ctrlpp::dare<Scalar, NX, NU>(sys_d.A, sys_d.B, Q, R).value().P;
    auto K = ctrlpp::lqr_gain<Scalar, NX, NU>(sys_d.A, sys_d.B, Q, R).value();

    std::printf("P_00,P_01,P_10,P_11,K_00,K_01\n");
    std::printf("%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n",
                P(0, 0), P(0, 1), P(1, 0), P(1, 1), K(0, 0), K(0, 1));
}
