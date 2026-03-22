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

    std::cout << "RLS Online Identification Example\n";
    std::cout << "True parameters: [" << true_params.transpose() << "]\n";
    std::cout << "Forgetting factor: " << lambda << "\n\n";

    for (int t = 1; t <= 200; ++t) {
        Eigen::Vector2d phi;
        phi << input(gen), input(gen);
        double y = true_params.dot(phi) + noise(gen);

        estimator.update(y, phi);

        if (t % 50 == 0) {
            auto theta = estimator.parameters();
            auto P = estimator.covariance();
            std::cout << "Step " << t << ":\n";
            std::cout << "  Parameters: [" << theta.transpose() << "]\n";
            std::cout << "  Covariance diagonal: [" << P.diagonal().transpose() << "]\n\n";
        }
    }

    auto final_params = estimator.parameters();
    std::cout << "Final parameters: [" << final_params.transpose() << "]\n";
    std::cout << "Parameter error:  ["
              << (final_params - true_params).transpose() << "]\n";

    return 0;
}
