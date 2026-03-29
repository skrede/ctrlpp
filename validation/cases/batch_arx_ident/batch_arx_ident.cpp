// batch_arx_ident.cpp -- Batch ARX identification matching batch_arx_ident.m
// Usage: ./batch_arx_ident > batch_arx_ident_cpp.csv

#include "ctrlpp/sysid/batch_arx.h"

#include <Eigen/Dense>

#include <cmath>
#include <cstdio>

int main()
{
    constexpr int n_samples = 200;

    Eigen::RowVectorXd u(n_samples);
    Eigen::RowVectorXd y(n_samples);
    u.setZero();
    y.setZero();

    for(int k = 0; k < n_samples; ++k)
        u(k) = std::sin(0.3 * (k + 1)) + 0.5 * std::cos(0.7 * (k + 1));

    // y(k) = 0.8*y(k-1) - 0.2*y(k-2) + 0.5*u(k-1) + 0.3*u(k-2)
    for(int k = 2; k < n_samples; ++k)
        y(k) = 0.8 * y(k - 1) - 0.2 * y(k - 2) + 0.5 * u(k - 1) + 0.3 * u(k - 2);

    auto result = ctrlpp::batch_arx<2, 2>(y, u);

    // Extract ARX parameters from observer canonical form state-space
    // A matrix first column contains a-coefficients (negated in some conventions)
    // We extract the raw theta vector by reading A and B
    auto& sys = result.system;

    // For ARX(2,2): the observer canonical form has
    // A = [[a1, 1], [a2, 0]], B = [[b1], [b2]], C = [1, 0]
    double a1 = sys.A(0, 0);
    double a2 = sys.A(1, 0);
    double b1 = sys.B(0, 0);
    double b2 = sys.B(1, 0);

    std::printf("a1,a2,b1,b2\n");
    std::printf("%.15e,%.15e,%.15e,%.15e\n", a1, a2, b1, b2);
}
