#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_TRAJECTORY_DEMO_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_TRAJECTORY_DEMO_H

#include "golden_reference.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/double_s_trajectory.h"

#include <cstdint>

namespace ctrlpp
{

inline const char* describe(trajectory_error e)
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
        case trajectory_error::boundary_velocity_rounding_ops_out_of_range:
            return "boundary velocity rounding count is out of range";
        case trajectory_error::unrepresentable_duration:
            return "duration is not representable";
        case trajectory_error::unreachable_boundary_velocity:
            return "boundary velocity is unreachable";
        case trajectory_error::acceleration_limit_infeasible:
            return "command needs more than the acceleration limit";
        case trajectory_error::duration_shorter_than_current:
            return "duration is shorter than the current one";
        case trajectory_error::unreachable_duration:
            return "duration is unreachable";
    }
    return "unknown";
}

// The instantiation the embedded compile witness builds for this family, so the family's allocation
// claims for construction and evaluation cover this demo.
template <class Scalar>
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

    // Dyadic limits that reach the velocity limit, so construction is closed form and exact at either
    // precision. The run's last sample, at 4 s of 4.125 s, lies in the final jerk segment.
    static typename profile_type::config command() { return {.q0 = Scalar{0}, .q1 = Scalar{6.75}, .v_max = Scalar{2}, .a_max = Scalar{4}, .j_max = Scalar{16}}; }

    Scalar step()
    {
        const Scalar t = static_cast<Scalar>(m_index) * static_cast<Scalar>(kDt);
        ++m_index;
        Scalar position{};
        Scalar velocity{};
        Scalar acceleration{};
        m_profile.evaluate(t, position, velocity, acceleration);
        return position;
    }

    const profile_type& profile() const { return m_profile; }

private:
    profile_type m_profile;
    std::int32_t m_index;

    explicit trajectory_demo(profile_type profile)
        : m_profile(profile)
        , m_index(0)
    {
    }
};

}

#endif
