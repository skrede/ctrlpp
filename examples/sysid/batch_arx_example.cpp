// Usage: ./batch_arx_example | gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '-' using 1:2 with lines title 'actual', '' using 1:3 with lines title 'predicted'"
// Redirect: ./batch_arx_example > output.csv
/// @file batch_arx_example.cpp
/// @brief Batch ARX system identification: identifies y(t) = 0.7*y(t-1) + 0.3*u(t-1)
///        and outputs actual vs predicted time series for gnuplot visualization.

#include "ctrlpp/sysid.h"

#include <Eigen/Dense>

#include <cstddef>
#include <iostream>
#include <random>

int main()
{
    // Batch ARX identification of y(t) = 0.7*y(t-1) + 0.3*u(t-1)
    constexpr std::size_t N = 500;

    Eigen::Matrix<double, 1, static_cast<int>(N)> Y;
    Eigen::Matrix<double, 1, static_cast<int>(N)> U;

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> u_dist(-1.0, 1.0);

    double y = 0.0;
    double u_prev = 0.0;
    for(std::size_t t = 0; t < N; ++t)
    {
        double u = u_dist(gen);
        double y_new = 0.7 * y + 0.3 * u_prev;
        Y(0, static_cast<int>(t)) = y_new;
        U(0, static_cast<int>(t)) = u;
        y = y_new;
        u_prev = u;
    }

    auto result = ctrlpp::batch_arx<1, 1>(Y, U);

    // Simulate the identified model to produce predicted output
    auto const& sys = result.system;
    Eigen::Vector<double, 1> x = Eigen::Vector<double, 1>::Zero();

    std::cout << "# step,actual,predicted\n";
    for(std::size_t t = 1; t < N; ++t)
    {
        Eigen::Vector<double, 1> u_t;
        u_t(0) = U(0, static_cast<int>(t - 1));
        x = sys.A * x + sys.B * u_t;
        Eigen::Vector<double, 1> y_hat = sys.C * x + sys.D * u_t;
        std::cout << t << ',' << Y(0, static_cast<int>(t)) << ',' << y_hat(0) << '\n';
    }

    // Print model info to stderr so it doesn't mix with CSV
    std::cerr << "Batch ARX: A=" << sys.A << " B=" << sys.B
              << " NRMSE=" << result.metrics.nrmse
              << " VAF=" << result.metrics.vaf << "%\n";

    return 0;
}
