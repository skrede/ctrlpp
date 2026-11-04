#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_PLANNER_3RD_DOUBLE_S_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_PLANNER_3RD_DOUBLE_S_H

#include "ctrlpp/trajectory/detail/planner_3rd_stop.h"
#include "ctrlpp/trajectory/detail/jerk_limited_ramp.h"

#include "ctrlpp/trajectory/double_s_trajectory.h"

#include <array>
#include <cmath>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar>
class planner_3rd_double_s : protected planner_3rd_stop<Scalar>
{
protected:
    using planner_3rd_stop<Scalar>::planner_3rd_stop;

    void append_double_s_phases(std::array<Scalar, 7> const& durations, Scalar sigma, Scalar peak)
    {
        auto const j_pos = sigma * this->m_j_max;
        auto const j_neg = -sigma * this->m_j_max;

        auto const ramp = this->m_phase_count;
        this->append_phase(durations[0], j_pos);
        append_hold_phase(durations[1], durations[2], peak, j_neg);
        this->append_phase(durations[2], j_neg);
        this->pin_group_end(ramp, peak, Scalar{0});
        auto const cruise = this->m_phase_count;
        this->append_phase(durations[3], Scalar{0});
        this->pin_group_end(cruise, peak, Scalar{0});
        this->m_target_anchor_first = this->m_phase_count - 1;
        auto const decel = this->m_phase_count;
        this->append_phase(durations[4], j_neg);
        append_hold_phase(durations[5], durations[6], Scalar{0}, j_pos);
        this->append_phase(durations[6], j_pos);
        this->pin_group_end(decel, Scalar{0}, Scalar{0});
    }

    // Pinned where the closing jerk phase's backward form from (v_end, 0) starts, so a sample in its second half
    // never passes v_end however short that phase is.
    void append_hold_phase(Scalar duration, Scalar closing_duration, Scalar v_end, Scalar j_close)
    {
        auto const hold = this->m_phase_count;
        this->append_phase(duration, Scalar{0});
        if(duration > Scalar{0})
        {
            auto const T_c = std::min(closing_duration, jerk_time_within(this->m_a_max, this->m_j_max));
            this->pin_group_end(hold, v_end + Scalar{0.5} * j_close * T_c * T_c, -j_close * T_c);
        }
    }

    // @cite biagiotti2009 -- Sec. 3.4.3
    void plan_rest_to_rest(Scalar q0)
    {
        auto const h_signed = this->m_target - q0;

        if(h_signed == Scalar{0})
        {
            return;
        }

        auto const sigma = (h_signed > Scalar{0}) ? Scalar{1} : Scalar{-1};
        auto const h = std::abs(h_signed);

        Scalar T_j1{}, T_a{}, T_v{}, T_d{}, T_j2{}, v_lim{};
        compute_double_s_durations(h, T_j1, T_a, T_v, T_d, T_j2, v_lim);

        append_double_s_phases({T_j1, T_a - Scalar{2} * T_j1, T_j1, T_v, T_j2, T_d - Scalar{2} * T_j2, T_j2}, sigma, sigma * v_lim);
    }

    // Neither v_max nor a_max reached.
    // @cite biagiotti2009 -- Sec. 3.4.3, p.90
    void assign_doubly_degenerate(Scalar h, Scalar& T_j1, Scalar& T_a, Scalar& T_d, Scalar& T_j2, Scalar& v_lim) const
    {
        auto const T_j_dd = doubly_degenerate_jerk_time(h, this->m_v_max, this->m_a_max, this->m_j_max);
        T_j1 = T_j_dd;
        T_j2 = T_j_dd;
        T_a = Scalar{2} * T_j_dd;
        T_d = Scalar{2} * T_j_dd;
        v_lim = this->m_j_max * T_j_dd * T_j_dd;
    }

    // c = v_max T_a - h >= 0 is the distance the cruise-free profile at v_max lacks. Each peak is formed at or below
    // v_max, by the same arithmetic as double_s_trajectory.
    // @cite biagiotti2009 -- Sec. 3.4.3, p.89-91
    void solve_no_cruise(Scalar h, Scalar c, bool a_max_reached, Scalar& T_j1, Scalar& T_a, Scalar& T_d, Scalar& T_j2, Scalar& v_lim) const
    {
        if(a_max_reached)
        {
            v_lim = no_cruise_peak(c, h, this->m_v_max, this->m_a_max, this->m_j_max);
            T_j1 = this->m_a_max / this->m_j_max;
            T_j2 = T_j1;
            T_a = T_j1 + v_lim / this->m_a_max;
            T_d = T_a;
            if(!(T_a < Scalar{2} * T_j1))
            {
                return;
            }
        }
        assign_doubly_degenerate(h, T_j1, T_a, T_d, T_j2, v_lim);
    }

    // @cite biagiotti2009 -- Sec. 3.4.3, p.88-91
    void compute_double_s_durations(Scalar h, Scalar& T_j1, Scalar& T_a, Scalar& T_v, Scalar& T_d, Scalar& T_j2, Scalar& v_lim) const
    {
        auto const v = this->m_v_max;
        auto const a = this->m_a_max;
        auto const j = this->m_j_max;

        bool const a_max_reached = (v * j >= a * a);

        auto const ramp = make_velocity_ramp(v, a, j);
        auto const c = v * ramp.T - h;

        if(c < Scalar{0})
        {
            T_j1 = ramp.T_j;
            T_j2 = ramp.T_j;
            T_a = ramp.T;
            T_d = ramp.T;
            T_v = -c / v;
            v_lim = this->m_v_max;
        }
        else
        {
            T_v = Scalar{0};
            solve_no_cruise(h, c, a_max_reached, T_j1, T_a, T_d, T_j2, v_lim);
        }
    }

    // Appends nothing when the shape does not exist for this command, so the caller can plan the move another way.
    // @cite biagiotti2009 -- Sec. 3.4.1, eq. (3.19)-(3.27), p.79-85
    bool append_incorporate_velocity(Scalar q0, Scalar v0)
    {
        auto const created = double_s_trajectory<Scalar>::create({
                .q0 = q0,
                .q1 = this->m_target,
                .v_max = this->m_v_max,
                .a_max = this->m_a_max,
                .j_max = this->m_j_max,
                .v0 = v0,
                .v1 = Scalar{0},
        });
        if(!created.has_value())
        {
            return false;
        }
        auto const& profile = *created;

        auto const durations = profile.phase_durations();
        auto const sigma = (this->m_target - q0 > Scalar{0}) ? Scalar{1} : Scalar{-1};
        append_double_s_phases(durations, sigma, sigma * profile.peak_velocity());
        return true;
    }
};

}

#endif
