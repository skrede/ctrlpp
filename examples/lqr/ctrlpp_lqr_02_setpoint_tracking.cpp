// Usage: ./ctrlpp_lqr_02_setpoint_tracking | gnuplot -p -e "plot '-' using 1:2 with lines"
// Redirect: ./ctrlpp_lqr_02_setpoint_tracking > output.csv

#include "ctrlpp/discretise_impl.h"
#include "ctrlpp/kalman.h"
#include "ctrlpp/lqr.h"
#include "ctrlpp/propagate.h"
#include "ctrlpp/state_space.h"

#include <iomanip>
#include <iostream>

int main()
{
    using Scalar = double;
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;
    constexpr std::size_t NY = 1;

    constexpr Scalar m = 1.0;
    constexpr Scalar k = 1.0;
    constexpr Scalar b = 0.5;
    constexpr Scalar dt = 0.05;
    constexpr Scalar duration = 20.0;

    ctrlpp::continuous_state_space<Scalar, NX, NU, NY> sys_c{};
    sys_c.A << 0.0, 1.0, -k / m, -b / m;
    sys_c.B << 0.0, 1.0 / m;
    sys_c.C << 1.0, 0.0;
    sys_c.D << 0.0;

    auto sys_d = ctrlpp::discretise(ctrlpp::zoh{}, sys_c, dt);

    Eigen::Matrix<Scalar, 2, 2> Q_lqr = Eigen::Matrix<Scalar, 2, 2>::Zero();
    Q_lqr(0, 0) = 10.0;
    Q_lqr(1, 1) = 1.0;
    Eigen::Matrix<Scalar, 1, 1> R_lqr;
    R_lqr << 1.0;

    auto K_opt = ctrlpp::lqr_gain<Scalar, NX, NU>(sys_d.A, sys_d.B, Q_lqr, R_lqr);
    ctrlpp::lqr<Scalar, NX, NU> controller(*K_opt);

    Eigen::Matrix<Scalar, 2, 2> I2 = Eigen::Matrix<Scalar, 2, 2>::Identity();
    Eigen::Matrix<Scalar, 1, 1> u_ff =
        sys_d.B.colPivHouseholderQr().solve((I2 - sys_d.A) * Eigen::Matrix<Scalar, 2, 1>{{2.0}, {0.0}});

    Eigen::Matrix<Scalar, 2, 2> Q_proc = Eigen::Matrix<Scalar, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<Scalar, 1, 1> R_meas;
    R_meas << 0.1;
    Eigen::Matrix<Scalar, 2, 2> P0 = Eigen::Matrix<Scalar, 2, 2>::Identity() * 100.0;

    Eigen::Matrix<Scalar, 2, 1> x0_est;
    x0_est << 5.0, -2.0;

    ctrlpp::kalman_filter<Scalar, NX, NU, NY> kf(sys_d, {.Q = Q_proc, .R = R_meas, .x0 = x0_est, .P0 = P0});

    Eigen::Matrix<Scalar, 2, 1> x_true = Eigen::Matrix<Scalar, 2, 1>::Zero();

    Eigen::Matrix<Scalar, 2, 1> x_ref;
    x_ref << 2.0, 0.0;

    std::cout << "time,x_true_0,x_true_1,x_est_0,x_est_1,x_ref_0,control\n";

    for (Scalar t = 0.0; t < duration; t += dt) {
        auto x_est = kf.state();
        Eigen::Matrix<Scalar, 1, 1> u = controller.compute(x_est - x_ref) + u_ff;

        Eigen::Matrix<Scalar, 1, 1> z = sys_d.C * x_true;

        std::cout << std::fixed << std::setprecision(4)
                  << t << ","
                  << x_true(0) << "," << x_true(1) << ","
                  << x_est(0) << "," << x_est(1) << ","
                  << x_ref(0) << ","
                  << u(0) << "\n";

        kf.predict(u);
        x_true = ctrlpp::propagate(sys_d, x_true, u);
        kf.update(z);
    }
}
