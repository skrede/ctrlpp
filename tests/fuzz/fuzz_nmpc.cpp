#include "ctrlpp/nmpc.h"
#include "ctrlpp/mpc/nlopt_solver.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

namespace {

struct double_integrator {
    static constexpr double dt = 0.1;

    auto operator()(const ctrlpp::Vector<double, 2>& x,
                    const ctrlpp::Vector<double, 1>& u) const
        -> ctrlpp::Vector<double, 2>
    {
        ctrlpp::Vector<double, 2> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    }
};

}

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // horizon_byte(1) + Q_diag(2*8=16) + R(8) + u_min(8) + u_max(8) +
    // x0(2*8=16) = 57 bytes
    if (size < 57)
        return 0;

    std::size_t offset = 0;

    // Horizon: 2-10
    int horizon = 2 + (data[offset] % 9);
    offset += 1;

    double buf[7];
    std::memcpy(buf, data + offset, 56);

    for (int i = 0; i < 7; ++i) {
        if (!std::isfinite(buf[i]))
            return 0;
    }

    // Q diagonal entries (positive via abs + epsilon, clamped)
    double q0 = std::min(std::abs(buf[0]) + 1e-6, 1e4);
    double q1 = std::min(std::abs(buf[1]) + 1e-6, 1e4);

    // R entry (positive, clamped)
    double r0 = std::min(std::abs(buf[2]) + 1e-6, 1e4);

    double u_min_val = buf[3];
    double u_max_val = buf[4];
    if (u_min_val >= u_max_val) {
        u_min_val = -10.0;
        u_max_val = 10.0;
    }

    ctrlpp::Vector<double, 2> x0;
    x0 << buf[5], buf[6];

    // Clamp x0 to prevent extreme values
    for (int i = 0; i < 2; ++i)
        x0(i) = std::clamp(x0(i), -1e3, 1e3);

    ctrlpp::nmpc_config<double, 2, 1> cfg;
    cfg.horizon = horizon;
    cfg.Q << q0, 0.0, 0.0, q1;
    cfg.R << r0;

    ctrlpp::Vector<double, 1> umin, umax;
    umin << u_min_val;
    umax << u_max_val;
    cfg.u_min = umin;
    cfg.u_max = umax;

    try {
        double_integrator dynamics;
        ctrlpp::nmpc<double, 2, 1, ctrlpp::nlopt_solver<double>, double_integrator>
            controller(dynamics, cfg);
        auto result = controller.solve(x0);

        if (result.has_value()) {
            const auto& u = result.value();
            for (int i = 0; i < 1; ++i) {
                if (!std::isfinite(u(i)))
                    __builtin_trap();
            }
        }
    } catch (...) {
        // NLopt may throw for degenerate inputs; that is acceptable
    }

    return 0;
}
