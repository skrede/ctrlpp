// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_lqr_04_pendulum_swingup' using 1:2 with lines title 'theta', '' using 1:3 with lines title 'theta_dot', '' using 1:4 with lines title 'torque', '' using 1:5 with lines title 'mode'"
// Redirect: ./ctrlpp_lqr_04_pendulum_swingup > output.csv

// This example self-verifies at the end of main and returns a nonzero exit code
// if the pendulum did not reach and hold upright.

#include "ctrlpp/control/lqr.h"
#include "ctrlpp/model/discretize.h"
#include "ctrlpp/model/state_space.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <numbers>
#include <iostream>
#include <algorithm>

// Classical two-mode pendulum swing-up: an energy-shaping law pumps the
// pendulum toward the upright energy level, then an LQR stabilizer catches and
// holds it near the top.
//
// Energy-shaping law after:
//   K. J. Astrom and K. Furuta, "Swinging up a pendulum by energy control",
//   Automatica 36(2), 2000, pp. 287-295.
//
// Model convention matches the NMPC swing-up example:
//   theta_ddot = +(g/l) sin(theta) - b*theta_dot + torque/(m l^2),
// so theta = 0 is the unstable upright equilibrium and theta = pi is the
// stable hanging one.
int main()
{
    using Scalar = double;
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;
    constexpr std::size_t NY = 2;

    constexpr Scalar g = 9.81;
    constexpr Scalar l = 0.5;
    constexpr Scalar m = 0.5;
    constexpr Scalar b = 0.1;
    constexpr Scalar dt = 0.05;
    constexpr Scalar duration = 12.0;

    // Torque authority is deliberately below the gravity torque m*g*l = 2.4525,
    // so the controller must pump energy over several swings instead of lifting
    // the pendulum directly.
    constexpr Scalar u_max = 2.0;

    // Energy at the upright equilibrium (theta = 0, theta_dot = 0).
    constexpr Scalar E_ref = m * g * l;

    // Energy-pump gain: sets how aggressively the shaping law drives the energy
    // error to zero. Positive by construction of the Astrom-Furuta law.
    constexpr Scalar k_e = 1.0;

    // Catch region around the upright: hand off to the LQR stabilizer once the
    // pendulum is both close to upright and slow enough to be captured.
    constexpr Scalar theta_catch = 0.35;  // rad, about 20 degrees from upright
    constexpr Scalar omega_catch = 3.0;   // rad/s

    // Nonlinear pendulum dynamics, forward-Euler integrated at dt.
    auto step = [&](const Eigen::Matrix<Scalar, NX, 1>& x, Scalar torque) -> Eigen::Matrix<Scalar, NX, 1>
    {
        const Scalar theta = x(0);
        const Scalar theta_dot = x(1);
        const Scalar theta_ddot = (g / l) * std::sin(theta) - b * theta_dot + torque / (m * l * l);
        return (Eigen::Matrix<Scalar, NX, 1>() << theta + theta_dot * dt, theta_dot + theta_ddot * dt).finished();
    };

    // Total mechanical energy with the sign consistent with the +(g/l) sin
    // dynamics: minimal (-m g l) hanging, maximal (+m g l) upright.
    auto energy = [&](const Eigen::Matrix<Scalar, NX, 1>& x) -> Scalar
    {
        return 0.5 * m * l * l * x(1) * x(1) + m * g * l * std::cos(x(0));
    };

    // Wrap an angle into (-pi, pi] so that "near upright" is measured against 0.
    auto wrap = [](Scalar angle) -> Scalar { return std::remainder(angle, 2.0 * std::numbers::pi_v<Scalar>); };

    // LQR stabilizer built from the discrete-time linearization about upright.
    // Linearizing theta_ddot about theta = 0 gives sin(theta) ~ theta, hence
    // A_c = [[0, 1], [g/l, -b]] and B_c = [[0], [1/(m l^2)]].
    ctrlpp::continuous_state_space<Scalar, NX, NU, NY> sys_c{};
    sys_c.A << 0.0, 1.0, g / l, -b;
    sys_c.B << 0.0, 1.0 / (m * l * l);
    sys_c.C.setIdentity();
    sys_c.D.setZero();

    auto sys_d = ctrlpp::discretize(ctrlpp::zoh{}, sys_c, dt);

    Eigen::Matrix<Scalar, NX, NX> Q_lqr = Eigen::Matrix<Scalar, NX, NX>::Zero();
    Q_lqr(0, 0) = 50.0;
    Q_lqr(1, 1) = 5.0;
    Eigen::Matrix<Scalar, NU, NU> R_lqr;
    R_lqr << 1.0;

    auto K_opt = ctrlpp::lqr_gain<Scalar, NX, NU>(sys_d.A, sys_d.B, Q_lqr, R_lqr);
    if(!K_opt)
    {
        std::cerr << "LQR gain synthesis failed\n";
        return EXIT_FAILURE;
    }
    ctrlpp::lqr<Scalar, NX, NU> stabilizer(*K_opt);

    // Start hanging at the stable equilibrium and swing up.
    Eigen::Matrix<Scalar, NX, 1> x(std::numbers::pi_v<Scalar>, 0.0);

    std::cout << "time,theta,theta_dot,torque,mode\n";

    for(Scalar t = 0.0; t < duration; t += dt)
    {
        const Scalar theta_wrapped = wrap(x(0));
        const Scalar theta_dot = x(1);

        Scalar torque = 0.0;
        int mode = 0;  // 0 = energy pump, 1 = LQR stabilizer

        if(std::abs(theta_wrapped) < theta_catch && std::abs(theta_dot) < omega_catch)
        {
            // Near upright: stabilize with LQR on the wrapped state.
            mode = 1;
            Eigen::Matrix<Scalar, NX, 1> x_lin(theta_wrapped, theta_dot);
            torque = stabilizer.compute(x_lin)(0);
        }
        else
        {
            // Energy-shaping law. The drive direction is sign(theta_dot cos theta);
            // ties (including the initial rest at the hanging point) resolve to +1
            // so the pendulum is kicked off the stable equilibrium.
            const Scalar drive = (theta_dot * std::cos(x(0)) >= 0.0) ? 1.0 : -1.0;
            torque = k_e * (energy(x) - E_ref) * drive;
        }

        torque = std::clamp(torque, -u_max, u_max);

        std::cout << std::fixed << std::setprecision(4) << t << "," << x(0) << "," << x(1) << "," << torque << "," << mode << "\n";

        x = step(x, torque);
    }

    // Self-check: the pendulum must have reached and be holding near upright.
    constexpr Scalar hold_theta = 0.1;  // rad
    constexpr Scalar hold_omega = 0.5;  // rad/s
    const Scalar theta_final = wrap(x(0));
    if(std::abs(theta_final) >= hold_theta || std::abs(x(1)) >= hold_omega)
    {
        std::cerr << "swing-up self-check FAILED: |theta|=" << std::abs(theta_final)
                  << " rad, |theta_dot|=" << std::abs(x(1)) << " rad/s\n";
        return EXIT_FAILURE;
    }

    std::cerr << "swing-up self-check passed: |theta|=" << std::abs(theta_final) << " rad, |theta_dot|=" << std::abs(x(1)) << " rad/s\n";

    return EXIT_SUCCESS;
}
