#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_TRAJECTORY_DEMO_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_TRAJECTORY_DEMO_H

#include "golden_reference.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/double_s_trajectory.h"

#include <cstdint>

namespace ctrlpp {

inline const char *describe(trajectory_error e)
{
    switch(e)
    {
        case trajectory_error::non_positive_velocity_limit:
            return "velocity limit is not positive";
        case trajectory_error::non_positive_acceleration_limit:
            return "acceleration limit is not positive";
        case trajectory_error::non_positive_jerk_limit:
            return "jerk limit is not positive";
        case trajectory_error::non_positive_duration:
            return "duration is not positive";
        case trajectory_error::non_finite_input:
            return "non-finite input";
        case trajectory_error::boundary_velocity_exceeds_limit:
            return "boundary velocity exceeds the limit";
        case trajectory_error::unrepresentable_duration:
            return "duration is not representable";
        case trajectory_error::unreachable_boundary_velocity:
            return "boundary velocity is unreachable";
        case trajectory_error::duration_shorter_than_current:
            return "duration is shorter than the current one";
        case trajectory_error::unreachable_duration:
            return "duration is unreachable";
    }
    return "unknown";
}

// double_s_trajectory from its fallible factory is what the embedded compile
// witness builds for this family, and the profile construction and evaluation
// rows of the allocation matrix make their claims about it.
template<class Scalar>
class trajectory_demo
{
public:
    using profile_type = double_s_trajectory<Scalar>;

    static ctrlpp::expected<trajectory_demo, trajectory_error> make()
    {
        const auto profile = profile_type::create(command());
        if(!profile.has_value())
            return ctrlpp::unexpected(profile.error());
        return trajectory_demo(*profile);
    }

    // Dyadic limits that reach the velocity limit, so construction is closed
    // form and exact at either precision. The run's last sample, at 4 s of
    // 4.125 s, lies in the final jerk segment.
    static typename profile_type::config command()
    {
        return {.q0 = Scalar{0}, .q1 = Scalar{6.75}, .v_max = Scalar{2}, .a_max = Scalar{4}, .j_max = Scalar{16}};
    }

    Scalar step()
    {
        const Scalar t = static_cast<Scalar>(index_) * static_cast<Scalar>(kDt);
        ++index_;
        Scalar position{};
        Scalar velocity{};
        Scalar acceleration{};
        profile_.evaluate(t, position, velocity, acceleration);
        return position;
    }

    const profile_type &profile() const
    {
        return profile_;
    }

private:
    profile_type profile_;
    std::int32_t index_;

    explicit trajectory_demo(profile_type profile)
            : profile_(profile)
            , index_(0)
    {
    }
};

}

#endif
