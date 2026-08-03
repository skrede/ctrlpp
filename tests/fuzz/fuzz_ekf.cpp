#include "ctrlpp/estimation/ekf.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>

namespace
{

// Bitwise identity of two matrices, which is what "left untouched" means and
// what `operator!=` cannot express. That operator is an elementwise IEEE
// comparison folded with `any()`, and IEEE says a NaN equals nothing, itself
// included -- so it reports a difference between a matrix and a byte-for-byte
// copy of itself the moment either carries one. A carried estimate that is
// already non-finite is precisely the state this check exists to inspect,
// because it is the only one that reaches the rejection path below at all.
template <typename Derived>
bool bitwise_equal(const Eigen::MatrixBase<Derived>& lhs, const Eigen::MatrixBase<Derived>& rhs)
{
    for(Eigen::Index i = 0; i < lhs.size(); ++i)
    {
        const auto a = lhs.derived().data()[i];
        const auto b = rhs.derived().data()[i];
        if(std::memcmp(&a, &b, sizeof(a)) != 0)
            return false;
    }
    return true;
}

}

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 10 doubles: x0(2), measurement(2), Q_diag(2), R_diag(2), A entries(2) = 80 bytes
    if(size < 80)
        return 0;

    double buf[10];
    std::memcpy(buf, data, 80);

    for(int i = 0; i < 10; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    Eigen::Matrix<double, 2, 1> x0;
    x0 << buf[0], buf[1];

    Eigen::Matrix<double, 2, 1> z;
    z << buf[2], buf[3];

    // Clamp covariance diagonals to positive
    double q0 = std::clamp(std::abs(buf[4]), 1e-10, 1e6);
    double q1 = std::clamp(std::abs(buf[5]), 1e-10, 1e6);
    double r0 = std::clamp(std::abs(buf[6]), 1e-10, 1e6);
    double r1 = std::clamp(std::abs(buf[7]), 1e-10, 1e6);

    // Clamp system matrix entries for stability
    double a00 = std::clamp(buf[8], -0.99, 0.99);
    double a11 = std::clamp(buf[9], -0.99, 0.99);

    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Zero();
    Q(0, 0) = q0;
    Q(1, 1) = q1;

    Eigen::Matrix<double, 2, 2> R = Eigen::Matrix<double, 2, 2>::Zero();
    R(0, 0) = r0;
    R(1, 1) = r1;

    // Nonlinear dynamics: x' = [a00*x0; a11*x1] (linear for simplicity, EKF handles it)
    auto dynamics = [a00, a11](const Eigen::Matrix<double, 2, 1>& x,
                               const Eigen::Matrix<double, 1, 1>& /*u*/) {
        Eigen::Matrix<double, 2, 1> xn;
        xn(0) = a00 * x(0);
        xn(1) = a11 * x(1);
        return xn;
    };

    // Identity measurement
    auto measurement = [](const Eigen::Matrix<double, 2, 1>& x) {
        return x;
    };

    ctrlpp::ekf_config<double, 2, 1, 2> cfg{.Q = Q, .R = R, .x0 = x0};
    // The harness refuses a non-finite raw input above and clamps every
    // configuration field, so the configuration validation cannot reject here. A
    // rejection would mean the validation refused a finite configuration, which
    // is a defect rather than a fuzz finding.
    auto created = ctrlpp::ekf<double, 2, 1, 2, decltype(dynamics), decltype(measurement)>::create(
        dynamics, measurement, cfg);
    if(!created)
        abort();
    auto& filter = *created;

    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();

    for(int step = 0; step < 10; ++step)
    {
        filter.predict(u);

        // The harness refuses a non-finite raw input above, so z is finite by
        // construction. A rejection can therefore only name the carried
        // estimate; a measurement rejection would mean the guard reported the
        // wrong cause, which is a defect rather than a fuzz finding. A rejected
        // step must also have left the estimate bitwise untouched, and the check
        // is a bitwise one: `predict` is infallible by contract and may poison
        // the covariance, so the state a rejection carries is routinely one an
        // IEEE comparison cannot even compare with itself.
        const Eigen::Matrix<double, 2, 1> x_before = filter.state();
        const Eigen::Matrix<double, 2, 2> P_before = filter.covariance();
        if(const auto stepped = filter.update(z); !stepped)
        {
            if(stepped.error() == ctrlpp::ekf_update_error::non_finite_measurement)
                abort();
            if(!bitwise_equal(filter.state(), x_before) || !bitwise_equal(filter.covariance(), P_before))
                abort();
            return 0;
        }

        auto const& x_est = filter.state();
        for(int i = 0; i < 2; ++i)
        {
            if(!std::isfinite(x_est(i)))
            return 0;
        }
    }

    return 0;
}
