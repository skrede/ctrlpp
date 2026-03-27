// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./rls_example' using 1:2 with lines title 'param 1', '' using 1:3 with lines title 'param 2'"
// Redirect: ./rls_example > output.csv
/// @file rls_example.cpp
/// @brief Online RLS identification: converges to true parameters [2.0, 0.5]
///        and outputs parameter estimates over time for gnuplot visualization.

#include "ctrlpp/sysid.h"

#include <Eigen/Dense>

#include <iostream>
#include <random>

int main()
{
    // Online RLS identification of y = 2.0*x1 + 0.5*x2
    constexpr std::size_t NP = 2;
    constexpr double lambda = 0.98;

    ctrlpp::rls_config<double, NP> cfg;
    cfg.lambda = lambda;
    ctrlpp::rls<double, NP> estimator(cfg);

    Eigen::Vector2d true_params;
    true_params << 2.0, 0.5;

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> input(-1.0, 1.0);
    std::normal_distribution<double> noise(0.0, 0.01);

    std::cout << "# step,param_1,param_2,cov_1,cov_2\n";

    for(int t = 1; t <= 200; ++t)
    {
        Eigen::Vector2d phi;
        phi << input(gen), input(gen);
        double y = true_params.dot(phi) + noise(gen);

        estimator.update(y, phi);

        auto theta = estimator.parameters();
        auto P = estimator.covariance();
        std::cout << t << ',' << theta(0) << ',' << theta(1) << ',' << P(0, 0) << ',' << P(1, 1) << '\n';
    }

    // Print final info to stderr
    auto final_params = estimator.parameters();
    std::cerr << "Final parameters: [" << final_params.transpose()
              << "] error: [" << (final_params - true_params).transpose() << "]\n";

    return 0;
}
