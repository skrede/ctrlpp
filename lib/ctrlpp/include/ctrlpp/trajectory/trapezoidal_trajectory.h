#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_TRAPEZOIDAL_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_TRAPEZOIDAL_TRAJECTORY_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
// Automatic Machines and Robots", 2009, Sec. 3.2, eq. (3.9)-(3.16), p.65-73

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/trapezoidal_types.h"
#include "ctrlpp/trajectory/trajectory_segment.h"

#include "ctrlpp/trajectory/detail/ramp_anchor.h"
#include "ctrlpp/trajectory/detail/boundary_velocity.h"
#include "ctrlpp/trajectory/detail/trapezoidal_profile.h"
#include "ctrlpp/trajectory/detail/trapezoidal_rescale.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <algorithm>

namespace ctrlpp
{

template <ctrlpp_floating_scalar Scalar>
class trapezoidal_trajectory : private detail::trapezoidal_rescale<Scalar>
{
    using base = detail::trapezoidal_rescale<Scalar>;

    struct frame_sample
    {
        Scalar q, dq, ddq;
    };

    struct unchecked_t
    {
        explicit unchecked_t() = default;
    };

public:
    using scalar_type = Scalar;

    struct config
    {
        Scalar q0, q1;
        Scalar v_max, a_max;
        Scalar v0{}, v1{};
        std::uint32_t boundary_velocity_rounding_ops{};
        acceleration_limit_policy acceleration_policy{acceleration_limit_policy::refuse};
    };

    // @cite biagiotti2009 -- Sec. 3.2.7, eq. (3.14)-(3.15), p.72
    static ctrlpp::expected<trapezoidal_trajectory, trajectory_error> create(config const& requested)
    {
        if(!std::isfinite(requested.q0) || !std::isfinite(requested.q1) || !std::isfinite(requested.v0) || !std::isfinite(requested.v1))
        {
            return ctrlpp::unexpected(trajectory_error::non_finite_input);
        }
        if(!std::isfinite(requested.v_max) || requested.v_max <= Scalar{0})
        {
            return ctrlpp::unexpected(trajectory_error::non_positive_velocity_limit);
        }
        if(!std::isfinite(requested.a_max) || requested.a_max <= Scalar{0})
        {
            return ctrlpp::unexpected(trajectory_error::non_positive_acceleration_limit);
        }
        auto const admitted = detail::admit_boundary_velocities(requested);
        if(!admitted.has_value())
        {
            return ctrlpp::unexpected(admitted.error());
        }
        if(!std::isfinite(detail::solve_trapezoidal_acceleration<Scalar>(*admitted)))
        {
            return ctrlpp::unexpected(trajectory_error::unreachable_boundary_velocity);
        }
        return checked(trapezoidal_trajectory{unchecked_t{}, *admitted}, *admitted);
    }

    // The clamp is defined only because `create` rejects a negative duration: std::clamp with its lower
    // bound above its upper bound is undefined behavior.
    // @cite biagiotti2009 -- Sec. 3.2.7, eq. (3.13a)-(3.13c), p.71
    trajectory_point<Scalar, 1> evaluate(Scalar t) const
    {
        auto const tc = std::clamp(t, Scalar{0}, this->m_duration);
        if(tc >= this->m_duration)
        {
            return {
                    .position = Vector<Scalar, 1>{this->m_q1},
                    .velocity = Vector<Scalar, 1>{this->m_sigma * this->m_v1},
                    .acceleration = Vector<Scalar, 1>{this->m_sigma * this->last_phase_acceleration()},
            };
        }
        // The stored total is a rounded sum, so the time back from it can exceed m_decel_duration at the switch;
        // each ramp's formula is used only strictly inside that ramp, timed from its own anchor.
        auto const dt_end = this->m_duration - tc;
        auto const s = (tc < this->m_accel_duration) ? acceleration_ramp(tc) : (dt_end >= this->m_decel_duration) ? cruise(tc) : deceleration_ramp(dt_end);
        return {
                .position = Vector<Scalar, 1>{this->m_q0 + this->m_sigma * s.q},
                .velocity = Vector<Scalar, 1>{this->m_sigma * s.dq},
                .acceleration = Vector<Scalar, 1>{this->m_sigma * s.ddq},
        };
    }

    Scalar duration() const { return this->m_duration; }

    // @cite biagiotti2009 -- Sec. 3.2.5, p.69
    bool is_triangular() const { return this->m_triangular; }

    Scalar peak_velocity() const { return this->m_sigma * this->m_v_v; }

    std::array<Scalar, 3> phase_durations() const { return {this->m_accel_duration, this->m_cruise_duration, this->m_decel_duration}; }

    trapezoidal_disposition<Scalar> const& disposition() const { return this->m_disposition; }

    // @cite biagiotti2009 -- Sec. 5.3, eq. (5.13)-(5.14) -- time scaling for synchronization
    ctrlpp::expected<void, trajectory_error> rescale_to(Scalar T_new)
    {
        auto const solved = this->solve_rescale(T_new);
        if(!solved.has_value())
        {
            return ctrlpp::unexpected(solved.error());
        }
        this->apply(*solved);
        return {};
    }

    // Runs the solve rescale_to runs and discards it, so a synchronized set can check every axis before it
    // commits any of them.
    ctrlpp::expected<void, trajectory_error> can_rescale_to(Scalar T_new) const
    {
        auto const solved = this->solve_rescale(T_new);
        if(!solved.has_value())
        {
            return ctrlpp::unexpected(solved.error());
        }
        return {};
    }

private:
    explicit trapezoidal_trajectory(unchecked_t, config const& cfg)
        : base{cfg}
    {
    }

    Scalar ramp_end_position() const { return this->m_v0 * this->m_accel_duration + Scalar{0.5} * this->m_a_a * this->m_accel_duration * this->m_accel_duration; }

    frame_sample acceleration_ramp(Scalar tc) const
    {
        auto const to_end = this->m_accel_duration - tc;
        return {detail::nearer_anchor_position(Scalar{0}, this->m_v0, ramp_end_position(), this->m_v_v, this->m_a_a, tc, to_end),
                detail::nearer_anchor_velocity(this->m_v0, this->m_v_v, this->m_a_a, tc, to_end), this->m_a_a};
    }

    frame_sample cruise(Scalar tc) const { return {ramp_end_position() + this->m_v_v * (tc - this->m_accel_duration), this->m_v_v, Scalar{0}}; }

    frame_sample deceleration_ramp(Scalar dt_end) const
    {
        auto const from_start = this->m_decel_duration - dt_end;
        auto const q_start = ramp_end_position() + this->m_v_v * this->m_cruise_duration;
        return {detail::nearer_anchor_position(q_start, this->m_v_v, std::abs(this->m_q1 - this->m_q0), this->m_v1, -this->m_a_d, from_start, dt_end),
                detail::nearer_anchor_velocity(this->m_v_v, this->m_v1, -this->m_a_d, from_start, dt_end), -this->m_a_d};
    }

    static ctrlpp::expected<trapezoidal_trajectory, trajectory_error> checked(trapezoidal_trajectory const& profile, config const& cfg)
    {
        for(auto const span : std::array{profile.m_accel_duration, profile.m_cruise_duration, profile.m_decel_duration, profile.m_duration})
        {
            if(span < Scalar{0})
            {
                return ctrlpp::unexpected(trajectory_error::unreachable_boundary_velocity);
            }
            if(!std::isfinite(span))
            {
                return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
            }
        }
        // A nonzero displacement or velocity change takes nonzero time at finite limits, so a total
        // that underflowed to zero describes no commanded motion. The predicate carries no threshold.
        if((cfg.q1 != cfg.q0 || cfg.v0 != cfg.v1) && !(profile.m_duration > Scalar{0}))
        {
            return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
        }
        if(cfg.acceleration_policy == acceleration_limit_policy::refuse && profile.m_a != cfg.a_max)
        {
            return ctrlpp::unexpected(trajectory_error::acceleration_limit_infeasible);
        }
        return profile;
    }
};

static_assert(trajectory_segment<trapezoidal_trajectory<double>, double, 1>);

}

#endif
