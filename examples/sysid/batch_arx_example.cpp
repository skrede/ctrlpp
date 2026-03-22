#include "ctrlpp/sysid.h"

#include <Eigen/Dense>

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
    for (std::size_t t = 0; t < N; ++t) {
        double u = u_dist(gen);
        double y_new = 0.7 * y + 0.3 * u_prev;
        Y(0, static_cast<int>(t)) = y_new;
        U(0, static_cast<int>(t)) = u;
        y = y_new;
        u_prev = u;
    }

    std::cout << "Batch ARX Identification Example\n";
    std::cout << "True system: y(t) = 0.7*y(t-1) + 0.3*u(t-1)\n\n";

    auto result = ctrlpp::batch_arx<1, 1>(Y, U);

    std::cout << "Identified state-space model:\n";
    std::cout << "  A = " << result.system.A << "\n";
    std::cout << "  B = " << result.system.B << "\n";
    std::cout << "  C = " << result.system.C << "\n";
    std::cout << "  D = " << result.system.D << "\n\n";

    std::cout << "Fit metrics:\n";
    std::cout << "  NRMSE: " << result.metrics.nrmse << "\n";
    std::cout << "  VAF:   " << result.metrics.vaf << "%\n\n";

    // Compare with online recursive ARX
    ctrlpp::recursive_arx<double, 1, 1> online_arx;
    for (std::size_t t = 0; t < N; ++t) {
        online_arx.update(Y(0, static_cast<int>(t)), U(0, static_cast<int>(t)));
    }

    auto online_ss = online_arx.to_state_space();
    std::cout << "Recursive ARX comparison:\n";
    std::cout << "  A = " << online_ss.A << "\n";
    std::cout << "  B = " << online_ss.B << "\n";

    return 0;
}
