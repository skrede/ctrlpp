#include "ctrlpp/estimation/kalman.h"
#include "ctrlpp/model/state_space.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // Ad(2x2=32) + Bd(2x1=16) + C(1x2=16) + D(1x1=8) + Q(2x2=32) + R(1x1=8)
    // + x0(2x1=16) + u(1x1=8) + z(1x1=8) = 144 bytes
    constexpr std::size_t needed = 144;
    if(size < needed)
        return 0;

    double buf[18];
    std::memcpy(buf, data, needed);

    for(int i = 0; i < 18; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    int idx = 0;

    // Clamp system matrix for stability (matching fuzz_ekf pattern)
    Eigen::Matrix<double, 2, 2> Ad = Eigen::Matrix<double, 2, 2>::Zero();
    Ad(0, 0) = std::clamp(buf[idx], -0.99, 0.99);
    Ad(1, 1) = std::clamp(buf[idx + 3], -0.99, 0.99);
    Ad(0, 1) = std::clamp(buf[idx + 1], -0.5, 0.5);
    Ad(1, 0) = std::clamp(buf[idx + 2], -0.5, 0.5);
    idx += 4;

    Eigen::Matrix<double, 2, 1> Bd;
    Bd << std::clamp(buf[idx], -1.0, 1.0), std::clamp(buf[idx + 1], -1.0, 1.0);
    idx += 2;

    Eigen::Matrix<double, 1, 2> C;
    C << std::clamp(buf[idx], -1.0, 1.0), std::clamp(buf[idx + 1], -1.0, 1.0);
    idx += 2;

    Eigen::Matrix<double, 1, 1> D;
    D << std::clamp(buf[idx], -1.0, 1.0);
    idx += 1;

    Eigen::Matrix<double, 2, 2> Q_raw;
    Q_raw << std::clamp(buf[idx], -1.0, 1.0), std::clamp(buf[idx + 1], -1.0, 1.0),
             std::clamp(buf[idx + 2], -1.0, 1.0), std::clamp(buf[idx + 3], -1.0, 1.0);
    idx += 4;

    // Make Q PSD: Q = Q_raw^T * Q_raw
    Eigen::Matrix<double, 2, 2> Q = Q_raw.transpose() * Q_raw;

    // Make R PD: R = r^2 + epsilon (clamp raw value to prevent overflow)
    double r_raw = std::clamp(buf[idx], -100.0, 100.0);
    double r_val = r_raw * r_raw + 1e-6;
    Eigen::Matrix<double, 1, 1> R;
    R << r_val;
    idx += 1;

    Eigen::Matrix<double, 2, 1> x0;
    x0 << std::clamp(buf[idx], -100.0, 100.0), std::clamp(buf[idx + 1], -100.0, 100.0);
    idx += 2;

    Eigen::Matrix<double, 1, 1> u;
    u << std::clamp(buf[idx], -100.0, 100.0);
    idx += 1;

    Eigen::Matrix<double, 1, 1> z;
    z << std::clamp(buf[idx], -100.0, 100.0);

    ctrlpp::discrete_state_space<double, 2, 1, 1> sys{Ad, Bd, C, D};
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity();

    ctrlpp::kalman_filter<double, 2, 1, 1> kf(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0});

    for(int step = 0; step < 10; ++step)
    {
        kf.predict(u);
        kf.update(z);

        const auto& x_est = kf.state();
        const auto& P_est = kf.covariance();

        for(int i = 0; i < 2; ++i)
        {
            if(!std::isfinite(x_est(i)))
                __builtin_trap();
            for(int j = 0; j < 2; ++j)
            {
                if(!std::isfinite(P_est(i, j)))
                    __builtin_trap();
            }
        }
    }

    return 0;
}
