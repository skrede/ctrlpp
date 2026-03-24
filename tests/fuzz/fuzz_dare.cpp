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

    Eigen::Matrix<double, 2, 2> A;
    A << buf[0], buf[1], buf[2], buf[3];

    Eigen::Matrix<double, 2, 1> B;
    B << buf[4], buf[5];

    Eigen::Matrix<double, 2, 2> Q_raw;
    Q_raw << buf[6], buf[7], buf[8], buf[9];

    // Make Q positive semi-definite: Q = Q_raw^T * Q_raw
    Eigen::Matrix<double, 2, 2> Q = Q_raw.transpose() * Q_raw;

    // Make R positive definite: R = R_raw^2 + epsilon
    double R_val = buf[10] * buf[10] + 1e-6;
    Eigen::Matrix<double, 1, 1> R;
    R << R_val;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);

    if(result.has_value())
    {
        const auto& P = result.value();
        for(int i = 0; i < 2; ++i)
        {
            for(int j = 0; j < 2; ++j)
            {
                if(!std::isfinite(P(i, j)))
                    __builtin_trap();
            }
        }
    }

    return 0;
}
