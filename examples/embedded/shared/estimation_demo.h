#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_ESTIMATION_DEMO_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_ESTIMATION_DEMO_H

#include "golden_reference.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/estimation/kalman.h"

#include <limits>
#include <cstdint>
#include <utility>

namespace ctrlpp {

inline const char *describe(filter_error e)
{
    switch(e)
    {
        case filter_error::degenerate_quaternion:
            return "initial quaternion is degenerate";
        case filter_error::non_positive_sigma_spread:
            return "sigma-point spread is not positive";
        case filter_error::non_positive_scaling_radicand:
            return "sigma-point scaling radicand is not positive";
        case filter_error::non_finite_process_noise:
            return "process noise is not finite";
        case filter_error::non_finite_measurement_noise:
            return "measurement noise is not finite";
        case filter_error::non_finite_initial_state:
            return "initial state is not finite";
        case filter_error::non_finite_initial_covariance:
            return "initial covariance is not finite";
    }
    return "unknown";
}

// kalman_filter at two states, one input and one output is the instantiation the
// embedded compile witness builds for this family, and the one the allocation
// matrix's estimation rows make their claim about.
template<class Scalar>
class estimation_demo
{
public:
    using filter_type     = kalman_filter<Scalar, kKalmanNx, 1, kKalmanNy>;
    using input_type      = typename filter_type::input_vector_t;
    using output_type     = typename filter_type::output_vector_t;
    using covariance_type = typename filter_type::cov_matrix_t;

    static ctrlpp::expected<estimation_demo, filter_error> make()
    {
        auto filter = filter_type::create(plant(), config());
        if(!filter.has_value())
            return ctrlpp::unexpected(filter.error());
        return estimation_demo(std::move(*filter));
    }

    // The compared quantity is the velocity variance, the largest entry the
    // covariance settles to. A refused update returns NaN, which no bound admits.
    Scalar step()
    {
        const Scalar position = static_cast<Scalar>(index_) * static_cast<Scalar>(kDt);
        ++index_;
        filter_.predict(input_type::Zero());
        if(!filter_.update(output_type{position}).has_value())
            return std::numeric_limits<Scalar>::quiet_NaN();
        return filter_.covariance()(1, 1);
    }

    const covariance_type &covariance() const
    {
        return filter_.covariance();
    }

    static typename filter_type::system_t plant()
    {
        const Scalar dt = static_cast<Scalar>(kDt);

        typename filter_type::system_t system{};
        system.A << Scalar{1}, dt, Scalar{0}, Scalar{1};
        system.B << Scalar{0.5} * dt * dt, dt;
        system.C << Scalar{1}, Scalar{0};
        system.D.setZero();
        return system;
    }

    // Dyadic entries, so the configuration rounds nothing at either precision.
    static kalman_config<Scalar, kKalmanNx, 1, kKalmanNy> config()
    {
        kalman_config<Scalar, kKalmanNx, 1, kKalmanNy> settings{};
        settings.Q << Scalar{1} / Scalar{1024}, Scalar{0}, Scalar{0}, Scalar{1} / Scalar{16};
        settings.R << Scalar{1};
        settings.x0.setZero();
        settings.P0 = covariance_type::Identity() / Scalar{16};
        return settings;
    }

private:
    filter_type filter_;
    std::int32_t index_;

    explicit estimation_demo(filter_type filter)
            : filter_(std::move(filter))
            , index_(0)
    {
    }
};

}

#endif
