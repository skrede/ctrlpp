#include "ctrlpp/control/l1.h"
#include "ctrlpp/control/l1_config.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 9 doubles: state x, reference r, predictor A, predictor B,
    //            gamma, bandwidth, sample_hz, theta_min, theta_max = 72 bytes
    if(size < 72)
        return 0;

    double buf[9];
    std::memcpy(buf, data, 72);

    for(int i = 0; i < 9; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    double x_val = std::clamp(buf[0], -1.0, 1.0);
    double r_val = std::clamp(buf[1], -1.0, 1.0);
    double pred_a = std::clamp(buf[2], -0.999, 0.999);
    double pred_b = std::clamp(buf[3], -10.0, 10.0);
    double gamma = std::clamp(std::abs(buf[4]), 1e-3, 1e6);
    double bandwidth = std::clamp(std::abs(buf[5]), 1.0, 1000.0);
    double sample_hz = std::clamp(std::abs(buf[6]), 100.0, 10000.0);
    double theta_min = -std::abs(buf[7]) - 1.0;
    double theta_max = std::abs(buf[8]) + 1.0;

    // Ensure bandwidth < Nyquist
    if(bandwidth >= sample_hz * 0.5)
        bandwidth = sample_hz * 0.4;

    ctrlpp::l1_config<double, 1, 1> cfg;
    cfg.predictor_model.A(0, 0) = pred_a;
    cfg.predictor_model.B(0, 0) = pred_b;
    cfg.predictor_model.C(0, 0) = 1.0;
    cfg.predictor_model.D(0, 0) = 0.0;
    cfg.gamma(0, 0) = gamma;
    cfg.theta_min(0) = theta_min;
    cfg.theta_max(0) = theta_max;

    ctrlpp::l1_controller<double, 1, 1> ctrl(cfg, bandwidth, sample_hz);

    Eigen::Matrix<double, 1, 1> x;
    Eigen::Matrix<double, 1, 1> r;
    r(0) = r_val;

    double plant_x = x_val;
    for(int step = 0; step < 20; ++step)
    {
        x(0) = plant_x;
        auto u = ctrl.evaluate(x, r);

        if(!std::isfinite(u(0)))
            __builtin_trap();

        // Simple plant: x[k+1] = 0.8 * x[k] + 0.5 * u[k]
        plant_x = 0.8 * plant_x + 0.5 * u(0);
    }

    return 0;
}
