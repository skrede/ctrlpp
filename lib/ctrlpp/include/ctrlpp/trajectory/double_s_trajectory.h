#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DOUBLE_S_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DOUBLE_S_TRAJECTORY_H

// @cite biagiotti2009 -- Sec. 3.4, eq. (3.17)-(3.33), Fig. 3.18, p.79-96

#include "ctrlpp/trajectory/trajectory_error.h"
#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/trajectory_segment.h"

#include "ctrlpp/trajectory/detail/double_s_rescale.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <algorithm>

namespace ctrlpp
{

template <ctrlpp_floating_scalar Scalar>
class double_s_trajectory : private detail::double_s_rescale<Scalar>
{
    using base = detail::double_s_rescale<Scalar>;

    struct unchecked_t
    {
        explicit unchecked_t() = default;
    };

public:
    using scalar_type = Scalar;

    struct config
    {
        Scalar q0, q1;
        Scalar v_max, a_max, j_max;
        Scalar v0{}, v1{};
        std::uint32_t boundary_velocity_rounding_ops{};
    };

    // @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2008, Sec. 3.4, eq. (3.17)-(3.18), p.80-81
    static ctrlpp::expected<double_s_trajectory, trajectory_error> create(config const& requested)
    {
        auto const limits = base::check_limits(requested);
        if(!limits.has_value())
        {
            return ctrlpp::unexpected(limits.error());
        }
        auto const admitted = detail::admit_boundary_velocities(requested);
        if(!admitted.has_value())
        {
            return ctrlpp::unexpected(admitted.error());
        }
        if(!base::is_reachable(*admitted))
        {
            return ctrlpp::unexpected(trajectory_error::unreachable_boundary_velocity);
        }
        return checked(double_s_trajectory{unchecked_t{}, *admitted}, *admitted);
    }

    // The clamp is defined only because `create` rejects a negative duration: std::clamp with its lower
    // bound above its upper bound is undefined behavior.
    // @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2008, Sec. 3.4.1, eq. (3.30a)-(3.30g), p.88-89
    trajectory_point<Scalar, 1> evaluate(Scalar t) const
    {
        if(this->m_duration <= Scalar{0})
        {
            return this->make_point(this->m_q0, Scalar{0}, Scalar{0});
        }
        auto const tc = std::clamp(t, Scalar{0}, this->m_duration);
        if(tc >= this->m_duration)
        {
            return this->make_point(this->m_q1, this->m_v1, Scalar{0});
        }
        if(tc <= Scalar{0})
        {
            return this->make_point(this->m_q0, this->m_v0, Scalar{0});
        }

        Scalar q{}, dq{}, ddq{};
        this->evaluate_positive_frame(tc, q, dq, ddq);

        // @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2008, Sec. 3.4.2, eq. (3.33), p.90
        auto const q_actual = this->m_q0 + this->m_sigma * (q - Scalar{0});
        auto const dq_actual = this->m_sigma * dq;
        auto const ddq_actual = this->m_sigma * ddq;
        return this->make_point(q_actual, dq_actual, ddq_actual);
    }

    // @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2008, Sec. 3.4.1, eq. (3.30a)-(3.30g), p.88-89
    void evaluate(Scalar t, Scalar& q, Scalar& dq, Scalar& ddq) const
    {
        if(this->m_duration <= Scalar{0})
        {
            this->store_point(q, dq, ddq, this->m_q0, Scalar{0}, Scalar{0});
            return;
        }
        auto const tc = std::clamp(t, Scalar{0}, this->m_duration);
        if(tc >= this->m_duration)
        {
            this->store_point(q, dq, ddq, this->m_q1, this->m_v1, Scalar{0});
            return;
        }
        if(tc <= Scalar{0})
        {
            this->store_point(q, dq, ddq, this->m_q0, this->m_v0, Scalar{0});
            return;
        }
        Scalar pq{}, pdq{}, pddq{};
        this->evaluate_positive_frame(tc, pq, pdq, pddq);
        // @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2008, Sec. 3.4.2, eq. (3.33), p.90
        this->store_point(q, dq, ddq, this->m_q0 + this->m_sigma * (pq - Scalar{0}), this->m_sigma * pdq, this->m_sigma * pddq);
    }

    Scalar duration() const { return this->m_duration; }

    bool is_degenerate() const { return this->m_degenerate; }

    Scalar peak_velocity() const { return this->m_v_lim; }

    // @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2008, Sec. 5.2, eq. (5.3)-(5.4), p.229-230
    ctrlpp::expected<void, trajectory_error> rescale_to(Scalar T_new)
    {
        auto const solved = this->solve_rescale(*this, T_new);
        if(!solved.has_value())
        {
            return ctrlpp::unexpected(solved.error());
        }
        *this = *solved;
        return {};
    }

    // Runs the solve rescale_to runs and discards it, so a synchronized set can check every axis before it
    // commits any of them.
    ctrlpp::expected<void, trajectory_error> can_rescale_to(Scalar T_new) const
    {
        auto const solved = this->solve_rescale(*this, T_new);
        if(!solved.has_value())
        {
            return ctrlpp::unexpected(solved.error());
        }
        return {};
    }

    // The floor absorbs the rounding of a constant-acceleration remainder that is algebraically zero at the
    // triangular boundary; a negative ramp duration is rejected by `create`, not hidden here.
    std::array<Scalar, 7> phase_durations() const
    {
        auto const const_accel = std::max(this->m_accel_duration - Scalar{2} * this->m_accel_jerk_duration, Scalar{0});
        auto const const_decel = std::max(this->m_decel_duration - Scalar{2} * this->m_decel_jerk_duration, Scalar{0});
        return {this->m_accel_jerk_duration, const_accel, this->m_accel_jerk_duration, this->m_cruise_duration, this->m_decel_jerk_duration, const_decel, this->m_decel_jerk_duration};
    }

private:
    explicit double_s_trajectory(unchecked_t, config const& cfg)
        : base{cfg}
    {
    }

    // The durations are checked as stored, not as phase_durations() reports them: that accessor floors
    // each constant-acceleration remainder at zero, and a negative ramp duration would vanish into the floor.
    static ctrlpp::expected<double_s_trajectory, trajectory_error> checked(double_s_trajectory const& profile, config const& cfg)
    {
        for(auto const stored : std::array<Scalar, 6>{
                    profile.m_accel_jerk_duration, profile.m_accel_duration, profile.m_cruise_duration, profile.m_decel_duration, profile.m_decel_jerk_duration, profile.m_duration})
        {
            if(stored < Scalar{0})
            {
                return ctrlpp::unexpected(trajectory_error::unreachable_boundary_velocity);
            }
            if(!std::isfinite(stored))
            {
                return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
            }
        }
        if(!(profile.m_v_lim <= cfg.v_max))
        {
            return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
        }
        // A nonzero displacement takes nonzero time at finite limits, so a total that underflowed to zero
        // describes no commanded motion. The predicate carries no threshold.
        if(cfg.q1 - cfg.q0 != Scalar{0} && !(profile.m_duration > Scalar{0}))
        {
            return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
        }
        return profile;
    }
};

static_assert(trajectory_segment<double_s_trajectory<double>, double, 1>);
static_assert(trajectory_segment<double_s_trajectory<float>, float, 1>);

}

#endif
