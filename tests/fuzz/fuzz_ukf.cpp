#include "ctrlpp/estimation/ukf.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 10 doubles: x0(2), measurement(2), Q_diag(2), R_diag(2), A entries(2) = 80 bytes
    if(size < 80)
        return 0;

    double buf[10];
    std::memcpy(buf, data, 80);

    for(int i = 0; i < 10; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    Eigen::Matrix<double, 2, 1> x0;
    x0 << buf[0], buf[1];

    Eigen::Matrix<double, 2, 1> z;
    z << buf[2], buf[3];

    double q0 = std::clamp(std::abs(buf[4]), 1e-10, 1e6);
    double q1 = std::clamp(std::abs(buf[5]), 1e-10, 1e6);
    double r0 = std::clamp(std::abs(buf[6]), 1e-10, 1e6);
    double r1 = std::clamp(std::abs(buf[7]), 1e-10, 1e6);

    double a00 = std::clamp(buf[8], -0.99, 0.99);
    double a11 = std::clamp(buf[9], -0.99, 0.99);

    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Zero();
    Q(0, 0) = q0;
    Q(1, 1) = q1;

    Eigen::Matrix<double, 2, 2> R = Eigen::Matrix<double, 2, 2>::Zero();
    R(0, 0) = r0;
    R(1, 1) = r1;

    auto dynamics = [a00, a11](const Eigen::Matrix<double, 2, 1>& x,
                               const Eigen::Matrix<double, 1, 1>& /*u*/) {
        Eigen::Matrix<double, 2, 1> xn;
        xn(0) = a00 * x(0);
        xn(1) = a11 * x(1);
        return xn;
    };

    auto measurement = [](const Eigen::Matrix<double, 2, 1>& x) {
        return x;
    };

    ctrlpp::ukf_config<double, 2, 1, 2> cfg{.Q = Q, .R = R, .x0 = x0};
    ctrlpp::ukf<double, 2, 1, 2, decltype(dynamics), decltype(measurement)> filter(
        dynamics, measurement, cfg);

    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();

    for(int step = 0; step < 10; ++step)
    {
        filter.predict(u);
        filter.update(z);

        auto const& x_est = filter.state();
        for(int i = 0; i < 2; ++i)
        {
            if(!std::isfinite(x_est(i)))
                __builtin_trap();
        }
    }

    return 0;
}
