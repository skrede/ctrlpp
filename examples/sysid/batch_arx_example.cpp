// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./batch_arx_example' using 1:2 with lines title 'actual', '' using 1:3 with lines title 'predicted'"
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

    // Identification is fallible: it rejects a record it cannot form a
    // regressor from, rather than returning a fit built on wrapped index
    // arithmetic. Handle the error branch explicitly.
    auto result = ctrlpp::batch_arx<1, 1>(Y, U);
    if(!result)
    {
        char const* reason = "unknown rejection";
        switch(result.error())
        {
        case ctrlpp::sysid_error::record_length_mismatch:
            reason = "the output and input records differ in length";
            break;
        case ctrlpp::sysid_error::record_not_single_row:
            reason = "a record is not a single row";
            break;
        case ctrlpp::sysid_error::too_few_samples:
            reason = "fewer samples than the model order requires";
            break;
        case ctrlpp::sysid_error::non_finite_sample:
            reason = "a sample is NaN or infinite";
            break;
        }
        std::cerr << "Batch ARX identification rejected the data: " << reason << '\n';
        return 1;
    }

    // Simulate the identified model to produce predicted output
    auto const& sys = result->system;
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
              << " NRMSE=" << result->metrics.nrmse
              << " VAF=" << result->metrics.vaf << "%\n";

    return 0;
}
