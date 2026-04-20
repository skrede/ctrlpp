// care_solution.cpp -- CARE solution matching care_solution.m
// Usage: ./care_solution > care_solution_cpp.csv

#include "ctrlpp/control/care.h"
#include "ctrlpp/control/lqr.h"

#include <cstdio>

int main()
{
    using Scalar = double;
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;

    Eigen::Matrix<Scalar, NX, NX> A;
    A << 0.0, 1.0, -1.0, -0.5;
    Eigen::Matrix<Scalar, NX, NU> B;
    B << 0.0, 1.0;

    Eigen::Matrix<Scalar, NX, NX> Q;
    Q << 10.0, 0.0, 0.0, 1.0;
    Eigen::Matrix<Scalar, NU, NU> R;
    R << 1.0;

    auto P = ctrlpp::care<Scalar, NX, NU>(A, B, Q, R).value().P;
    auto K = ctrlpp::lqr_gain_continuous<Scalar, NX, NU>(A, B, Q, R).value();

    std::printf("P_00,P_01,P_10,P_11,K_00,K_01\n");
    std::printf("%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n",
                P(0, 0), P(0, 1), P(1, 0), P(1, 1), K(0, 0), K(0, 1));
}
