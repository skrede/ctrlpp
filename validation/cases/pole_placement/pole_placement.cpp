// pole_placement.cpp -- Pole placement matching pole_placement.m
// Usage: ./pole_placement > pole_placement_cpp.csv

#include "ctrlpp/model/discretise.h"
#include "ctrlpp/model/state_space.h"
#include "ctrlpp/control/place.h"

#include <Eigen/Eigenvalues>

#include <array>
#include <complex>
#include <cstdio>
#include <iostream>

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

    std::array<std::complex<Scalar>, 2> desired = {
        std::complex<Scalar>{0.5, 0.1},
        std::complex<Scalar>{0.5, -0.1}
    };

    auto K_result = ctrlpp::place<Scalar, NX, NU>(sys_d.A, sys_d.B, desired);
    if(!K_result.has_value())
    {
        std::cerr << "pole placement declined the validation plant\n";
        return 1;
    }
    auto K = *K_result;

    Eigen::Matrix<Scalar, 2, 2> Acl = sys_d.A - sys_d.B * K;
    Eigen::EigenSolver<Eigen::Matrix<Scalar, 2, 2>> es(Acl);
    auto cl_eig = es.eigenvalues();

    // Sort eigenvalues by imaginary part (positive first) to match Octave ordering
    int i0 = 0, i1 = 1;
    if(cl_eig(0).imag() < cl_eig(1).imag())
    {
        i0 = 1;
        i1 = 0;
    }

    std::printf("K_00,K_01,cl_eig0_re,cl_eig0_im,cl_eig1_re,cl_eig1_im\n");
    std::printf("%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n",
                K(0, 0), K(0, 1),
                cl_eig(i0).real(), cl_eig(i0).imag(),
                cl_eig(i1).real(), cl_eig(i1).imag());
}
