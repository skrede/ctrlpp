// Usage: ./ctrlpp_lqr_observer_console | gnuplot -p -e "plot '-' using 1:3 with lines"
// Redirect: ./ctrlpp_lqr_observer_console > output.csv

#include "ctrlpp/ctrlpp.h"
#include "ctrlpp/ctrlpp_eigen.h"

#include <iomanip>
#include <iostream>

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;
    constexpr std::size_t NY = 1;
    constexpr double dt = 0.1;

    using Scalar = double;
    using Mat2 = Eigen::Matrix<Scalar, 2, 2>;
    using Vec2 = Eigen::Matrix<Scalar, 2, 1>;
    using Mat1 = Eigen::Matrix<Scalar, 1, 1>;
    using MatC = Eigen::Matrix<Scalar, 1, 2>;
    using MatB = Eigen::Matrix<Scalar, 2, 1>;
    using MatD = Eigen::Matrix<Scalar, 1, 1>;
    using System = ctrlpp::DiscreteStateSpace<Scalar, NX, NU, NY, ctrlpp::EigenLinalgPolicy>;

    Mat2 A;
    A << 1.0, dt,
         0.0, 1.0;
    MatB B;
    B << 0.5 * dt * dt,
         dt;
    MatC C;
    C << 1.0, 0.0;
    MatD D = MatD::Zero();
    System sys{A, B, C, D};

    Mat2 Q = Mat2::Zero();
    Q(0, 0) = 10.0;
    Q(1, 1) = 1.0;
    Mat1 R;
    R << 1.0;

    auto K_opt = ctrlpp::lqr_gain<Scalar, NX, NU>(A, B, Q, R);
    if (!K_opt) return 1;
    auto lqr = ctrlpp::Lqr<Scalar, NX, NU>(*K_opt);

    Mat2 Q_proc = Mat2::Identity() * 0.01;
    Mat1 R_meas;
    R_meas << 0.1;
    Mat2 P0 = Mat2::Identity();
    Vec2 x0_est = Vec2::Zero();

    auto kf = ctrlpp::KalmanFilter<Scalar, NX, NU, NY>(sys, Q_proc, R_meas, x0_est, P0);

    Vec2 x_true;
    x_true << 1.0, 0.0;
    Mat1 u_ctrl = Mat1::Zero();
    constexpr double duration = 10.0;

    std::cout << "time,setpoint,measurement,control,state_0,state_1,est_0,est_1" << std::endl;

    for (double t = 0.0; t < duration; t += dt) {
        kf.predict(u_ctrl);
        x_true = (A * x_true + B * u_ctrl).eval();
        Mat1 z = (C * x_true).eval();
        kf.update(z);
        u_ctrl = lqr.compute(kf.state());

        std::cout << std::fixed << std::setprecision(4)
                  << t << "," << 0.0 << "," << x_true(0) << "," << u_ctrl(0, 0)
                  << "," << x_true(0) << "," << x_true(1)
                  << "," << kf.state()(0) << "," << kf.state()(1) << std::endl;
    }
}
