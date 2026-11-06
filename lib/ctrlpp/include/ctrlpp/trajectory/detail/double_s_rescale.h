#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_DOUBLE_S_RESCALE_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_DOUBLE_S_RESCALE_H

#include "ctrlpp/trajectory/detail/double_s_profile.h"

#include "ctrlpp/trajectory/trajectory_error.h"

#include "ctrlpp/expected.h"

#include <cmath>
#include <algorithm>

namespace ctrlpp::detail
{

// The profile type is a template parameter because a rescale rebuilds the profile through its own `create`;
// `self` is the profile this base belongs to.
template <typename Scalar>
class double_s_rescale : protected double_s_profile<Scalar>
{
protected:
    using double_s_profile<Scalar>::double_s_profile;

    // With both boundary velocities at rest the duration is proportional to the reciprocal of the scale, so
    // one rebuild settles it; otherwise the scale is found by bisection, terminating on bracket exhaustion.
    // @cite biagiotti2009 -- Sec. 5.3 -- time scaling for multi-axis synchronization
    template <typename Profile>
    ctrlpp::expected<Profile, trajectory_error> solve_rescale(Profile const& self, Scalar T_new) const
    {
        // A synchronized set passes the slowest axis a bit-exact copy of its own duration, so equality is
        // tested first and that axis never reaches the shortening rejection.
        if(T_new == this->m_duration)
        {
            return self;
        }
        if(!std::isfinite(T_new) || T_new <= Scalar{0})
        {
            return ctrlpp::unexpected(trajectory_error::non_positive_duration);
        }
        if(T_new < this->m_duration)
        {
            return ctrlpp::unexpected(trajectory_error::duration_shorter_than_current);
        }
        if(!(this->m_duration > Scalar{0}))
        {
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }
        auto const cfg = command<typename Profile::config>();
        if(this->m_v0 == Scalar{0} && this->m_v1 == Scalar{0})
        {
            return rebuild_scaled<Profile>(cfg, this->m_duration / T_new);
        }
        return bisect_rescale<Profile>(cfg, T_new);
    }

    template <typename Config>
    Config command() const
    {
        return Config{.q0 = this->m_q0, .q1 = this->m_q1, .v_max = this->m_v_max, .a_max = this->m_a_max, .j_max = this->m_j_max, .v0 = this->m_v0, .v1 = this->m_v1};
    }

    template <typename Profile>
    ctrlpp::expected<Profile, trajectory_error> bisect_rescale(typename Profile::config const& cfg, Scalar T_new) const
    {
        auto const lambda_min = std::max(std::abs(this->m_v0), std::abs(this->m_v1)) / this->m_v_max;
        if(!(lambda_min < Scalar{1}))
        {
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }
        // lo and hi end as neighbors with the crossing between them; when lo is inadmissible or falls short,
        // the crossing lies in the range no admissible profile covers.
        auto const solved = rebuild_scaled<Profile>(cfg, bisect_scale<Profile>(cfg, lambda_min, T_new));
        if(!solved.has_value() || !(solved->duration() >= T_new))
        {
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }
        return solved;
    }

    template <typename Profile>
    static Scalar bisect_scale(typename Profile::config const& cfg, Scalar lo, Scalar T_new)
    {
        auto hi = Scalar{1};
        for(;;)
        {
            auto const mid = lo + (hi - lo) / Scalar{2};
            if(!(mid > lo) || !(mid < hi))
            {
                break;
            }
            if(reaches_duration<Profile>(cfg, mid, T_new))
            {
                lo = mid;
            }
            else
            {
                hi = mid;
            }
        }
        return lo;
    }

    // A rebuilt profile is admissible only with a positive duration and a peak at or below its own scaled
    // velocity limit; the boundary velocities are the command, not a limit, and stay unscaled.
    // @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2008, Sec. 5.2, eq. (5.3)-(5.4), p.229-230
    template <typename Profile>
    static ctrlpp::expected<Profile, trajectory_error> rebuild_scaled(typename Profile::config const& cfg, Scalar lambda)
    {
        auto scaled = cfg;
        scaled.v_max = cfg.v_max * lambda;
        scaled.a_max = cfg.a_max * lambda * lambda;
        scaled.j_max = cfg.j_max * lambda * lambda * lambda;

        // A single return of the named result builds it in the caller's slot; returning a separate error
        // object beside it grew the Cortex-M7 rescale frames.
        auto rebuilt = Profile::create(scaled);
        if(!rebuilt.has_value() || !(rebuilt->duration() > Scalar{0}) || !(rebuilt->peak_velocity() <= scaled.v_max))
        {
            rebuilt = ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }
        return rebuilt;
    }

    // A scale with no admissible profile lies below the admissible range and answers true, and inside it the
    // duration falls as the scale grows, so the answer is monotone and bracket halving is valid.
    template <typename Profile>
    static bool reaches_duration(typename Profile::config const& cfg, Scalar lambda, Scalar T_new)
    {
        auto const rebuilt = rebuild_scaled<Profile>(cfg, lambda);
        return !rebuilt.has_value() || rebuilt->duration() >= T_new;
    }
};

}

#endif
