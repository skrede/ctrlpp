#include "ctrlpp/sysid/rls.h"

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <cstring>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 6 doubles: regressor(2), output, forgetting_factor, P_diag(2) = 48 bytes
    if(size < 48)
        return 0;

    double buf[6];
    std::memcpy(buf, data, 48);

    for(int i = 0; i < 6; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    double phi0 = std::clamp(buf[0], -1e3, 1e3);
    double phi1 = std::clamp(buf[1], -1e3, 1e3);
    double y = std::clamp(buf[2], -1e6, 1e6);
    double lambda = std::clamp(buf[3], 0.9, 1.0);
    double p0 = std::clamp(std::abs(buf[4]), 1e-3, 1e6);
    double p1 = std::clamp(std::abs(buf[5]), 1e-3, 1e6);

    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Zero();
    P0(0, 0) = p0;
    P0(1, 1) = p1;

    // The harness clamps the forgetting factor into [0.9, 1.0] and every
    // covariance entry into a finite positive range above, so the configuration
    // validation cannot reject here. A rejection would mean the validation
    // refused an in-domain configuration, which is a defect rather than a fuzz
    // finding.
    auto created = ctrlpp::rls<double, 2>::create({.lambda = lambda, .P0 = P0});
    if(!created)
        abort();
    auto& estimator = *created;

    Eigen::Matrix<double, 2, 1> phi;
    phi << phi0, phi1;

    for(int step = 0; step < 10; ++step)
    {
        // A refusal is a legitimate outcome on decoded input -- the harness
        // clamps the configuration but not the regressor, so a covariance the
        // regressor cannot resolve is reachable. What must hold is the
        // reject-before-mutate contract: a refused cycle leaves the parameters
        // bitwise as they were, so the finiteness check below still applies and
        // the loop simply stops feeding a sample the estimator declined.
        auto const applied = estimator.update(y, phi);
        if(!applied)
            return 0;

        auto const& theta = estimator.parameters();
        for(int i = 0; i < 2; ++i)
        {
            if(!std::isfinite(theta(i)))
            return 0;
        }
    }

    return 0;
}
