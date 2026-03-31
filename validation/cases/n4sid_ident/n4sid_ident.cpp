// n4sid_ident.cpp -- N4SID identification matching n4sid_ident.m
// Usage: ./n4sid_ident > n4sid_ident_cpp.csv

#include "ctrlpp/sysid/n4sid.h"

#include <Eigen/Dense>

#include <cmath>
#include <cstdio>

int main()
{
    constexpr int n_samples = 500;

    // Generate data from same system as Octave script
    Eigen::Matrix2d A_true;
    A_true << 0.8, 0.1, -0.2, 0.9;
    Eigen::Vector2d B_true{0.5, 0.3};
    Eigen::RowVector2d C_true{1.0, 0.0};

    Eigen::RowVectorXd u(n_samples);
    Eigen::RowVectorXd y(n_samples);
    Eigen::Vector2d x = Eigen::Vector2d::Zero();

    for(int k = 0; k < n_samples; ++k)
    {
        u(k) = std::sin(0.3 * (k + 1)) + 0.5 * std::cos(0.7 * (k + 1)) + 0.3 * std::sin(1.1 * (k + 1));
        y(k) = (C_true * x)(0);
        x = A_true * x + B_true * u(k);
    }

    auto result = ctrlpp::n4sid<2>(y, u);

    auto& sys = result.system;

    // Simulate identified model
    Eigen::Vector2d x_id = Eigen::Vector2d::Zero();

    std::printf("step,y_actual,y_predicted\n");

    for(int k = 0; k < n_samples; ++k)
    {
        double y_pred = (sys.C * x_id + sys.D * Eigen::Matrix<double, 1, 1>::Constant(u(k)))(0);
        std::printf("%.15e,%.15e,%.15e\n",
                    static_cast<double>(k + 1), y(k), y_pred);
        Eigen::Matrix<double, 1, 1> u_vec;
        u_vec << u(k);
        x_id = (sys.A * x_id + sys.B * u_vec).eval();
    }
}
