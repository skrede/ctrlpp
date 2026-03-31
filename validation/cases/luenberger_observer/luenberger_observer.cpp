// luenberger_observer.cpp -- Luenberger observer matching luenberger_observer.m
// Usage: ./luenberger_observer > luenberger_observer_cpp.csv

#include "ctrlpp/estimation/luenberger.h"
#include "ctrlpp/model/discretise.h"
#include "ctrlpp/model/propagate.h"
#include "ctrlpp/model/state_space.h"
#include "ctrlpp/control/place.h"

#include <array>
#include <cmath>
#include <complex>
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

    constexpr Scalar dt = 0.05;
    constexpr Scalar duration = 5.0;

    auto sys_d = ctrlpp::discretise(ctrlpp::zoh{}, sys_c, dt);

    std::array<std::complex<Scalar>, 2> desired_obs = {
        std::complex<Scalar>{0.3, 0.0},
        std::complex<Scalar>{0.2, 0.0}
    };
    auto L = ctrlpp::place_observer<Scalar, NX, NY>(sys_d.A, sys_d.C, desired_obs).value();

    Eigen::Matrix<Scalar, 2, 1> x0_est = Eigen::Matrix<Scalar, 2, 1>::Zero();
    ctrlpp::luenberger_observer<Scalar, NX, NU, NY> obs(sys_d, L, x0_est);

    Eigen::Matrix<Scalar, 2, 1> x_true;
    x_true << 1.0, 0.0;

    std::printf("time,x_true_0,x_true_1,x_est_0,x_est_1\n");

    for(Scalar t = 0.0; t < duration - dt / 2; t += dt)
    {
        Eigen::Matrix<Scalar, 1, 1> u;
        u << 0.5 * std::sin(t);

        auto x_est = obs.state();
        std::printf("%.15e,%.15e,%.15e,%.15e,%.15e\n",
                    t, x_true(0), x_true(1), x_est(0), x_est(1));

        obs.predict(u);
        x_true = ctrlpp::propagate(sys_d, x_true, u);

        Eigen::Matrix<Scalar, 1, 1> z_new = sys_d.C * x_true;
        obs.update(z_new);
    }
}
