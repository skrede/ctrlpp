#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/state_space.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // horizon_byte(1) + Q_diag(2*8=16) + R(8) + u_min(8) + u_max(8) +
    // x0(2*8=16) + x_ref(2*8=16) = 73 bytes
    if (size < 73)
        return 0;

    std::size_t offset = 0;

    // Horizon: 2-20
    int horizon = 2 + (data[offset] % 19);
    offset += 1;

    double buf[9];
    std::memcpy(buf, data + offset, 72);

    for (int i = 0; i < 9; ++i) {
        if (!std::isfinite(buf[i]))
            return 0;
    }

    // Q diagonal entries (positive via abs + epsilon)
    double q0 = std::abs(buf[0]) + 1e-6;
    double q1 = std::abs(buf[1]) + 1e-6;

    // R entry (positive)
    double r0 = std::abs(buf[2]) + 1e-6;

    // Clamp Q, R to reasonable range to avoid numerical blowup
    q0 = std::min(q0, 1e6);
    q1 = std::min(q1, 1e6);
    r0 = std::min(r0, 1e6);

    double u_min_val = buf[3];
    double u_max_val = buf[4];
    if (u_min_val >= u_max_val) {
        u_min_val = -1.0;
        u_max_val = 1.0;
    }

    Eigen::Matrix<double, 2, 1> x0;
    x0 << buf[5], buf[6];

    Eigen::Matrix<double, 2, 1> x_ref;
    x_ref << buf[7], buf[8];

    // Fixed double-integrator system (discrete, dt=0.1)
    constexpr double dt = 0.1;
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, dt, 0.0, 1.0;

    Eigen::Matrix<double, 2, 1> B;
    B << 0.5 * dt * dt, dt;

    Eigen::Matrix<double, 2, 2> C = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 2, 1> D = Eigen::Matrix<double, 2, 1>::Zero();

    ctrlpp::discrete_state_space<double, 2, 1, 2> sys{A, B, C, D};

    ctrlpp::mpc_config<double, 2, 1> cfg;
    cfg.horizon = horizon;
    cfg.Q << q0, 0.0, 0.0, q1;
    cfg.R << r0;

    Eigen::Matrix<double, 1, 1> umin, umax;
    umin << u_min_val;
    umax << u_max_val;
    cfg.u_min = umin;
    cfg.u_max = umax;

    try {
        ctrlpp::mpc<double, 2, 1, ctrlpp::osqp_solver> controller(sys, cfg);
        auto result = controller.solve(x0, x_ref);

        if (result.has_value()) {
            const auto& u = result.value();
            for (int i = 0; i < 1; ++i) {
                if (!std::isfinite(u(i)))
                    __builtin_trap();
            }
        }
    } catch (...) {
        // OSQP setup or solve may throw for degenerate inputs; that is acceptable
    }

    return 0;
}
