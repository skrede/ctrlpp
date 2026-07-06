#include "ctrlpp/control/pid.h"
#include "ctrlpp/control/pid_policies.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 8 doubles: kp, ki, kd, output_min, output_max, dt, setpoint, measurement = 64 bytes
    if(size < 64)
        return 0;

    double buf[8];
    std::memcpy(buf, data, 64);

    for(int i = 0; i < 8; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    double kp = std::clamp(buf[0], -100.0, 100.0);
    double ki = std::clamp(buf[1], -100.0, 100.0);
    double kd = std::clamp(buf[2], -100.0, 100.0);
    double output_min = std::clamp(buf[3], -1e3, 1e3);
    double output_max = std::clamp(buf[4], -1e3, 1e3);
    double dt = buf[5];
    double setpoint = std::clamp(buf[6], -100.0, 100.0);
    double measurement = std::clamp(buf[7], -100.0, 100.0);

    // Ensure dt > 0
    if(dt < 1e-6)
        dt = 1e-6;

    // Ensure output_min < output_max
    if(output_min >= output_max)
    {
        double tmp = output_min;
        output_min = output_max - 1.0;
        output_max = tmp + 1.0;
        if(output_min >= output_max)
            return 0;
    }

    using PidType = ctrlpp::pid<double, 1, ctrlpp::anti_windup<ctrlpp::clamping>>;
    using ConfigType = PidType::config_type;

    ConfigType cfg;
    cfg.kp[0] = kp;
    cfg.ki[0] = ki;
    cfg.kd[0] = kd;
    cfg.output_min[0] = output_min;
    cfg.output_max[0] = output_max;

    PidType pid(cfg);

    Eigen::Matrix<double, 1, 1> sp;
    sp << setpoint;
    Eigen::Matrix<double, 1, 1> meas;
    meas << measurement;

    for(int step = 0; step < 20; ++step)
    {
        auto u = pid.compute(sp, meas, dt);
        // NaN from integral overflow with extreme gains is expected behavior.
        // The harness only checks that the library does not segfault.
        if(!std::isfinite(u(0)))
            return 0;
    }

    return 0;
}
