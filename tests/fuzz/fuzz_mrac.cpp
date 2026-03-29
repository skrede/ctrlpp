#include "ctrlpp/control/mrac.h"
#include "ctrlpp/control/mrac_config.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 8 doubles: state x, reference r, gamma, ref_model_pole, sign_b, dt,
    //            theta_x_init, theta_r_init = 64 bytes
    if(size < 64)
        return 0;

    double buf[8];
    std::memcpy(buf, data, 64);

    for(int i = 0; i < 8; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    double x_val = std::clamp(buf[0], -1e3, 1e3);
    double r_val = std::clamp(buf[1], -1e3, 1e3);
    double gamma = std::clamp(std::abs(buf[2]), 1e-3, 1e6);
    double a_m = std::clamp(buf[3], -0.999, 0.999);
    double sign_b = (buf[4] >= 0.0) ? 1.0 : -1.0;
    double dt = std::clamp(std::abs(buf[5]), 1e-6, 1.0);
    double theta_x0 = std::clamp(buf[6], -1e3, 1e3);
    double theta_r0 = std::clamp(buf[7], -1e3, 1e3);

    // SISO MRAC with dead-zone robustification
    using MracType = ctrlpp::mrac_controller<double, 1, 1, ctrlpp::dead_zone>;
    using Config = ctrlpp::mrac_config<double, 1, 1, ctrlpp::dead_zone>;

    // Reference model: x_m[k+1] = a_m * x_m[k] + (1 - a_m) * r[k]
    Config cfg;
    cfg.reference_model.A(0, 0) = a_m;
    cfg.reference_model.B(0, 0) = 1.0 - std::abs(a_m);
    cfg.reference_model.C(0, 0) = 1.0;
    cfg.reference_model.D(0, 0) = 0.0;
    cfg.gamma_x(0, 0) = gamma;
    cfg.gamma_r(0, 0) = gamma;
    cfg.sign_b(0, 0) = sign_b;
    cfg.theta_x_0(0, 0) = theta_x0;
    cfg.theta_r_0(0, 0) = theta_r0;
    cfg.robustification.threshold = 0.01;

    MracType mrac(cfg);

    Eigen::Matrix<double, 1, 1> x;
    Eigen::Matrix<double, 1, 1> r;
    r(0) = r_val;

    double plant_x = x_val;
    for(int step = 0; step < 20; ++step)
    {
        x(0) = plant_x;
        auto u = mrac.evaluate(x, r);

        if(!std::isfinite(u(0)))
            __builtin_trap();

        // Simple plant: x[k+1] = 0.8 * x[k] + 0.5 * u[k]
        plant_x = 0.8 * plant_x + 0.5 * u(0);
    }

    return 0;
}
