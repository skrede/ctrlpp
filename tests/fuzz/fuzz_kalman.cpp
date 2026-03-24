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
    if (size < needed)
        return 0;

    double buf[18];
    std::memcpy(buf, data, needed);

    for (int i = 0; i < 18; ++i) {
        if (!std::isfinite(buf[i]))
            return 0;
    }

    int idx = 0;

    Eigen::Matrix<double, 2, 2> Ad;
    Ad << buf[idx], buf[idx + 1], buf[idx + 2], buf[idx + 3];
    idx += 4;

    Eigen::Matrix<double, 2, 1> Bd;
    Bd << buf[idx], buf[idx + 1];
    idx += 2;

    Eigen::Matrix<double, 1, 2> C;
    C << buf[idx], buf[idx + 1];
    idx += 2;

    Eigen::Matrix<double, 1, 1> D;
    D << buf[idx];
    idx += 1;

    Eigen::Matrix<double, 2, 2> Q_raw;
    Q_raw << buf[idx], buf[idx + 1], buf[idx + 2], buf[idx + 3];
    idx += 4;

    // Make Q PSD: Q = Q_raw^T * Q_raw
    Eigen::Matrix<double, 2, 2> Q = Q_raw.transpose() * Q_raw;

    // Make R PD: R = r^2 + epsilon
    double r_val = buf[idx] * buf[idx] + 1e-6;
    Eigen::Matrix<double, 1, 1> R;
    R << r_val;
    idx += 1;

    Eigen::Matrix<double, 2, 1> x0;
    x0 << buf[idx], buf[idx + 1];
    idx += 2;

    Eigen::Matrix<double, 1, 1> u;
    u << buf[idx];
    idx += 1;

    Eigen::Matrix<double, 1, 1> z;
    z << buf[idx];

    ctrlpp::discrete_state_space<double, 2, 1, 1> sys{Ad, Bd, C, D};
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity();

    ctrlpp::kalman_filter<double, 2, 1, 1> kf(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0});

    for (int step = 0; step < 10; ++step) {
        kf.predict(u);
        kf.update(z);

        const auto& x_est = kf.state();
        const auto& P_est = kf.covariance();

        for (int i = 0; i < 2; ++i) {
            if (!std::isfinite(x_est(i)))
                __builtin_trap();
            for (int j = 0; j < 2; ++j) {
                if (!std::isfinite(P_est(i, j)))
                    __builtin_trap();
            }
        }
    }

    return 0;
}
