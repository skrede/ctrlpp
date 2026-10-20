#include "ctrlpp/sysid/moesp.h"

#include <Eigen/Dense>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <algorithm>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 20 doubles: 10 input/output pairs (u_i, y_i) = 160 bytes
    constexpr int n_samples = 10;
    constexpr std::size_t needed = n_samples * 2 * sizeof(double);
    if(size < needed)
        return 0;

    double buf[n_samples * 2];
    std::memcpy(buf, data, needed);

    bool const finite = std::all_of(buf, buf + n_samples * 2, [](double v) { return std::isfinite(v); });

    Eigen::Matrix<double, 1, Eigen::Dynamic> U(1, n_samples);
    Eigen::Matrix<double, 1, Eigen::Dynamic> Y(1, n_samples);

    for(int i = 0; i < n_samples; ++i)
    {
        U(0, i) = finite ? std::clamp(buf[i * 2], -1e3, 1e3) : buf[i * 2];
        Y(0, i) = finite ? std::clamp(buf[i * 2 + 1], -1e3, 1e3) : buf[i * 2 + 1];
    }

    // Ten samples with three block rows pass every shape precondition at order
    // two: (3 - 1) * 1 >= 2 and 10 - 6 + 1 >= 3 + 2. The only refusals left are
    // a non-finite sample and an identification that produced no model.
    auto const result = ctrlpp::moesp<2>(Y, U, 3);

    if(!finite)
    {
        if(result.has_value() || result.error() != ctrlpp::sysid_error::non_finite_sample)
            abort();
        return 0;
    }

    if(!result)
    {
        if(result.error() != ctrlpp::sysid_error::model_not_identified)
            abort();
        return 0;
    }

    auto const& sys = result->system;
    if(!sys.A.allFinite() || !sys.B.allFinite() || !sys.C.allFinite() || !sys.D.allFinite())
        abort();
    if(!std::isfinite(result->condition_number))
        abort();

    return 0;
}
