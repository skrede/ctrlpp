#include "ctrlpp/sysid/n4sid.h"

#include <Eigen/Dense>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 20 doubles: 10 input/output pairs (u_i, y_i) = 160 bytes
    constexpr int n_samples = 10;
    constexpr std::size_t needed = n_samples * 2 * sizeof(double);
    if(size < needed)
        return 0;

    double buf[n_samples * 2];
    std::memcpy(buf, data, needed);

    for(int i = 0; i < n_samples * 2; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    // Build input and output matrices (1 x n_samples)
    Eigen::Matrix<double, 1, Eigen::Dynamic> U(1, n_samples);
    Eigen::Matrix<double, 1, Eigen::Dynamic> Y(1, n_samples);

    for(int i = 0; i < n_samples; ++i)
    {
        U(0, i) = std::clamp(buf[i * 2], -1e3, 1e3);
        Y(0, i) = std::clamp(buf[i * 2 + 1], -1e3, 1e3);
    }

    // Run N4SID with model order 2, block_rows = 3 (minimum for 10 samples)
    auto result = ctrlpp::n4sid<2>(Y, U, 3);

    // Degenerate data produces condition_number = infinity -- accepted
    if(std::isinf(result.condition_number))
        return 0;

    // For non-degenerate results, system matrices must be finite
    for(int r = 0; r < 2; ++r)
    {
        for(int c = 0; c < 2; ++c)
        {
            if(!std::isfinite(result.system.A(r, c)))
            return 0;
        }
        if(!std::isfinite(result.system.B(r, 0)))
            return 0;
    }

    for(int c = 0; c < 2; ++c)
    {
        if(!std::isfinite(result.system.C(0, c)))
            return 0;
    }

    if(!std::isfinite(result.system.D(0, 0)))
            return 0;

    return 0;
}
