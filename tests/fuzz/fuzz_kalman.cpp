#include "ctrlpp/estimation/kalman.h"

#include "ctrlpp/model/state_space.h"

#include "ctrlpp/detail/covariance_ops.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <optional>

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

    // Every field is finite, but Q is a Gram matrix: symmetric positive
    // semidefinite in exact arithmetic, while its rounded product can leave the
    // triangles apart by rounding drift or put an eigenvalue below the
    // definiteness floor. The factory must admit a matrix exactly when its own
    // copy symmetrizes within drift and the projected copy is positive
    // semidefinite, and otherwise name the first failing check in declaration
    // order, asymmetry ahead of definiteness within a matrix; anything else is
    // a defect rather than a fuzz finding.
    auto verdict = [](auto M, ctrlpp::filter_error asymmetric, ctrlpp::filter_error indefinite) -> std::optional<ctrlpp::filter_error> {
        if(!ctrlpp::detail::symmetrize_within_drift(M))
            return asymmetric;
        if(!ctrlpp::detail::psd_decomposition(M).psd)
            return indefinite;
        return std::nullopt;
    };
    const auto q_fault = verdict(Q, ctrlpp::filter_error::asymmetric_process_noise, ctrlpp::filter_error::indefinite_process_noise);
    const auto r_fault = verdict(R, ctrlpp::filter_error::asymmetric_measurement_noise, ctrlpp::filter_error::indefinite_measurement_noise);
    const auto p_fault = verdict(P0, ctrlpp::filter_error::asymmetric_initial_covariance, ctrlpp::filter_error::indefinite_initial_covariance);
    const auto first = q_fault ? q_fault : r_fault ? r_fault : p_fault;
    auto created = ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0});
    if(created.has_value() == first.has_value())
        abort();
    if(!created)
    {
        if(created.error() != *first)
            abort();
        return 0;
    }
    auto& kf = *created;

    for(int step = 0; step < 10; ++step)
    {
        kf.predict(u);

        // The harness refuses a non-finite raw input above, so z is finite by
        // construction. A rejection can therefore only name the carried
        // estimate; a measurement rejection would mean the guard reported the
        // wrong cause, which is a defect rather than a fuzz finding. A rejected
        // step must also have left the estimate bitwise untouched.
        const Eigen::Matrix<double, 2, 1> x_before = kf.state();
        const Eigen::Matrix<double, 2, 2> P_before = kf.covariance();
        if(const auto stepped = kf.update(z); !stepped)
        {
            if(stepped.error() == ctrlpp::kalman_update_error::non_finite_measurement)
                abort();
            if(kf.state() != x_before || kf.covariance() != P_before)
                abort();
            return 0;
        }

        const auto& x_est = kf.state();
        const auto& P_est = kf.covariance();

        for(int i = 0; i < 2; ++i)
        {
            if(!std::isfinite(x_est(i)))
            return 0;
            for(int j = 0; j < 2; ++j)
            {
                if(!std::isfinite(P_est(i, j)))
            return 0;
            }
        }
    }

    return 0;
}
