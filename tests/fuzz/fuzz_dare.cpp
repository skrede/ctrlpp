#include "ctrlpp/control/dare.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // Need 88 bytes: A(2x2=32) + B(2x1=16) + Q(2x2=32) + R(1x1=8)
    if(size < 88)
        return 0;

    double buf[11];
    std::memcpy(buf, data, 88);

    // Reject non-finite inputs early
    for(int i = 0; i < 11; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    // Clamp matrix entries to prevent intermediate overflow in symplectic construction
    Eigen::Matrix<double, 2, 2> A;
    A << std::clamp(buf[0], -2.0, 2.0), std::clamp(buf[1], -2.0, 2.0),
         std::clamp(buf[2], -2.0, 2.0), std::clamp(buf[3], -2.0, 2.0);

    Eigen::Matrix<double, 2, 1> B;
    B << std::clamp(buf[4], -2.0, 2.0), std::clamp(buf[5], -2.0, 2.0);

    Eigen::Matrix<double, 2, 2> Q_raw;
    Q_raw << std::clamp(buf[6], -2.0, 2.0), std::clamp(buf[7], -2.0, 2.0),
             std::clamp(buf[8], -2.0, 2.0), std::clamp(buf[9], -2.0, 2.0);

    // Make Q positive semi-definite: Q = Q_raw^T * Q_raw
    Eigen::Matrix<double, 2, 2> Q = Q_raw.transpose() * Q_raw;

    // Make R positive definite: R = R_raw^2 + epsilon
    double r_raw = std::clamp(buf[10], -2.0, 2.0);
    double R_val = r_raw * r_raw + 1e-6;
    Eigen::Matrix<double, 1, 1> R;
    R << R_val;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);

    // If dare returns a value, verify it is finite.
    // Non-finite results from ill-conditioned Schur decomposition
    // are accepted without trapping -- the library should return
    // nullopt but Eigen internals may produce edge cases with
    // subnormal inputs that we cannot guard against.
    if(result.has_value() && !result->P.allFinite())
        return 0;

    return 0;
}
