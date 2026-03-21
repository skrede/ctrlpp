// Usage: ./ctrlpp_lqr_03_observer_comparison | gnuplot -p -e "plot '-' using 1:6 with lines"
// Redirect: ./ctrlpp_lqr_03_observer_comparison > output.csv

#include "ctrlpp/discretise_impl.h"
#include "ctrlpp/kalman.h"
#include "ctrlpp/luenberger.h"
#include "ctrlpp/lqr.h"
#include "ctrlpp/propagate.h"
#include "ctrlpp/state_space.h"

#include <iomanip>
#include <iostream>

int main()
{
    using Scalar = double;
    constexpr std::size_t NX = 4;
    constexpr std::size_t NU = 1;
    constexpr std::size_t NY = 2;

    constexpr Scalar M_cart = 1.0;
    constexpr Scalar m_pend = 0.1;
    constexpr Scalar l = 0.5;
    constexpr Scalar g = 9.81;
    constexpr Scalar dt = 0.05;
    constexpr Scalar duration = 20.0;

    constexpr Scalar denom = M_cart + m_pend;

    ctrlpp::continuous_state_space<Scalar, NX, NU, NY> sys_c{};
    sys_c.A << 0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, -m_pend * g / denom, 0.0,
               0.0, 0.0, 0.0, 1.0,
               0.0, 0.0, (denom * g) / (denom * l), 0.0;
    sys_c.B << 0.0, 1.0 / denom, 0.0, -1.0 / (denom * l);
    sys_c.C << 1.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 1.0, 0.0;
    sys_c.D.setZero();

    auto sys_d = ctrlpp::discretise(ctrlpp::zoh{}, sys_c, dt);

    Eigen::Matrix<Scalar, 4, 4> Q_lqr = Eigen::Matrix<Scalar, 4, 4>::Zero();
    Q_lqr(0, 0) = 10.0;
    Q_lqr(1, 1) = 1.0;
    Q_lqr(2, 2) = 100.0;
    Q_lqr(3, 3) = 10.0;
    Eigen::Matrix<Scalar, 1, 1> R_lqr;
    R_lqr << 1.0;

    auto K_opt = ctrlpp::lqr_gain<Scalar, NX, NU>(sys_d.A, sys_d.B, Q_lqr, R_lqr);
    ctrlpp::lqr<Scalar, NX, NU> controller(*K_opt);

    Eigen::Matrix<Scalar, 4, 4> Q_proc = Eigen::Matrix<Scalar, 4, 4>::Identity() * 0.01;
    Eigen::Matrix<Scalar, 2, 2> R_meas = Eigen::Matrix<Scalar, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<Scalar, 4, 4> P0 = Eigen::Matrix<Scalar, 4, 4>::Identity();
    Eigen::Matrix<Scalar, 4, 1> x0_est = Eigen::Matrix<Scalar, 4, 1>::Zero();

    ctrlpp::kalman_filter<Scalar, NX, NU, NY> kf(sys_d, Q_proc, R_meas, x0_est, P0);

    Eigen::Matrix<Scalar, 4, 4> Q_obs = Eigen::Matrix<Scalar, 4, 4>::Identity() * 100.0;
    Eigen::Matrix<Scalar, 2, 2> R_obs = Eigen::Matrix<Scalar, 2, 2>::Identity();
    auto L_dual = ctrlpp::lqr_gain<Scalar, NX, NY>(
        sys_d.A.transpose(), sys_d.C.transpose(), Q_obs, R_obs);
    Eigen::Matrix<Scalar, 4, 2> L = L_dual->transpose();

    ctrlpp::luenberger_observer<Scalar, NX, NU, NY> luen(sys_d, L, x0_est);

    Eigen::Matrix<Scalar, 4, 1> x_true;
    x_true << 0.1, 0.0, 0.05, 0.0;

    std::cout << "time,"
              << "x_true_0,x_true_1,x_true_2,x_true_3,"
              << "kalman_0,kalman_1,kalman_2,kalman_3,"
              << "luenberger_0,luenberger_1,luenberger_2,luenberger_3,"
              << "control\n";

    for (Scalar t = 0.0; t < duration; t += dt) {
        auto x_kf = kf.state();
        auto x_lu = luen.state();
        auto u = controller.compute(x_kf);

        Eigen::Matrix<Scalar, 2, 1> z = sys_d.C * x_true;

        std::cout << std::fixed << std::setprecision(4)
                  << t << ","
                  << x_true(0) << "," << x_true(1) << ","
                  << x_true(2) << "," << x_true(3) << ","
                  << x_kf(0) << "," << x_kf(1) << ","
                  << x_kf(2) << "," << x_kf(3) << ","
                  << x_lu(0) << "," << x_lu(1) << ","
                  << x_lu(2) << "," << x_lu(3) << ","
                  << u(0) << "\n";

        kf.predict(u);
        luen.predict(u);
        x_true = ctrlpp::propagate(sys_d, x_true, u);
        kf.update(z);
        luen.update(z);
    }
}
