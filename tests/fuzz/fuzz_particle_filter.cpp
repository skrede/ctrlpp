#include "ctrlpp/estimation/particle_filter.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <random>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 8 doubles: measurement(2), Q_diag(2), R_diag(2), x0(2) = 64 bytes
    if(size < 64)
        return 0;

    double buf[8];
    std::memcpy(buf, data, 64);

    for(int i = 0; i < 8; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    Eigen::Matrix<double, 2, 1> z;
    z << std::clamp(buf[0], -1e3, 1e3), std::clamp(buf[1], -1e3, 1e3);

    double q0 = std::clamp(std::abs(buf[2]), 1e-10, 1e4);
    double q1 = std::clamp(std::abs(buf[3]), 1e-10, 1e4);
    double r0 = std::clamp(std::abs(buf[4]), 1e-10, 1e4);
    double r1 = std::clamp(std::abs(buf[5]), 1e-10, 1e4);

    Eigen::Matrix<double, 2, 1> x0;
    x0 << std::clamp(buf[6], -1e3, 1e3), std::clamp(buf[7], -1e3, 1e3);

    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Zero();
    Q(0, 0) = q0;
    Q(1, 1) = q1;

    Eigen::Matrix<double, 2, 2> R = Eigen::Matrix<double, 2, 2>::Zero();
    R(0, 0) = r0;
    R(1, 1) = r1;

    auto dynamics = [](const Eigen::Matrix<double, 2, 1>& x,
                       const Eigen::Matrix<double, 1, 1>& /*u*/) {
        Eigen::Matrix<double, 2, 1> xn;
        xn(0) = 0.9 * x(0);
        xn(1) = 0.9 * x(1);
        return xn;
    };

    auto measurement = [](const Eigen::Matrix<double, 2, 1>& x) {
        return x;
    };

    ctrlpp::pf_config<double, 2, 1, 2> cfg{.Q = Q, .R = R, .x0 = x0};

    // Use deterministic RNG seed for reproducibility
    std::mt19937_64 rng(42);

    ctrlpp::particle_filter<double, 2, 1, 2, 50, decltype(dynamics), decltype(measurement)>
        filter(dynamics, measurement, cfg, rng);

    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();

    for(int step = 0; step < 10; ++step)
    {
        filter.predict(u);
        filter.update(z);

        auto x_est = filter.state();
        for(int i = 0; i < 2; ++i)
        {
            if(!std::isfinite(x_est(i)))
                __builtin_trap();
        }
    }

    return 0;
}
