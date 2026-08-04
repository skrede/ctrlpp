// analysis_poles.cpp -- System analysis matching analysis_poles.m
// Usage: ./analysis_poles > analysis_poles_cpp.csv

#include "ctrlpp/model/analysis.h"
#include "ctrlpp/model/discretize.h"
#include "ctrlpp/model/state_space.h"

#include <algorithm>
#include <cstdio>

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

    auto sys_d = ctrlpp::discretize(ctrlpp::zoh{}, sys_c, 0.05);

    auto pc = ctrlpp::poles(sys_c);
    auto pd = ctrlpp::poles(sys_d);

    // Sort by imaginary part descending
    std::sort(pc.begin(), pc.end(), [](auto a, auto b) { return a.imag() > b.imag(); });
    std::sort(pd.begin(), pd.end(), [](auto a, auto b) { return a.imag() > b.imag(); });

    bool ctrb = ctrlpp::is_controllable<Scalar, NX, NU>(sys_d.A, sys_d.B);
    bool obsv = ctrlpp::is_observable<Scalar, NX, NY>(sys_d.A, sys_d.C);
    bool stable_c = ctrlpp::is_stable(sys_c);
    bool stable_d = ctrlpp::is_stable(sys_d);

    std::printf("pc0_re,pc0_im,pc1_re,pc1_im,pd0_re,pd0_im,pd1_re,pd1_im,ctrb_rank,obsv_rank,stable_c,stable_d\n");
    std::printf("%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n",
                pc[0].real(), pc[0].imag(), pc[1].real(), pc[1].imag(),
                pd[0].real(), pd[0].imag(), pd[1].real(), pd[1].imag(),
                ctrb ? 2.0 : 1.0, obsv ? 2.0 : 1.0,
                stable_c ? 1.0 : 0.0, stable_d ? 1.0 : 0.0);
}
