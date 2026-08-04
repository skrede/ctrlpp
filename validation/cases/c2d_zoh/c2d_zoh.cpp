// c2d_zoh.cpp -- ZOH discretisation matching c2d_zoh.m
// Usage: ./c2d_zoh > c2d_zoh_cpp.csv

#include "ctrlpp/model/discretize.h"
#include "ctrlpp/model/state_space.h"

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

    auto sys_d = ctrlpp::discretize(ctrlpp::zoh{}, sys_c, 0.05);

    std::printf("Ad_00,Ad_01,Ad_10,Ad_11,Bd_00,Bd_10\n");
    std::printf("%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n",
                sys_d.A(0, 0), sys_d.A(0, 1), sys_d.A(1, 0), sys_d.A(1, 1),
                sys_d.B(0, 0), sys_d.B(1, 0));
}
