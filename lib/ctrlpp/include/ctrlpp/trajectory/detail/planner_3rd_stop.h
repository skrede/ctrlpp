#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_PLANNER_3RD_STOP_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_PLANNER_3RD_STOP_H

#include "ctrlpp/trajectory/detail/jerk_limited_ramp.h"
#include "ctrlpp/trajectory/detail/planner_3rd_state.h"

#include <cmath>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar>
class planner_3rd_stop : protected planner_3rd_state<Scalar>
{
protected:
    struct stop_result
    {
        Scalar dist; // signed
        Scalar T_j;
        Scalar brake_duration;
        bool a_max_reached;
    };

    using planner_3rd_state<Scalar>::planner_3rd_state;

    stop_result compute_stop(Scalar v0) const
    {
        auto const abs_v = std::abs(v0);
        auto const sign_v = (v0 >= Scalar{0}) ? Scalar{1} : Scalar{-1};

        Scalar T_j{};
        Scalar T_total{};
        Scalar dist{};
        bool a_reached{};

        if(abs_v > this->m_a_max * this->m_a_max / this->m_j_max)
        {
            T_j = jerk_time_within(this->m_a_max, this->m_j_max);
            auto const T_const = abs_v / this->m_a_max - T_j;
            T_total = T_j + T_const + T_j;
            dist = sign_v * compute_decel_distance(abs_v, T_j, T_const);
            a_reached = true;
        }
        else
        {
            T_j = std::min(std::sqrt(abs_v / this->m_j_max), jerk_time_within(this->m_a_max, this->m_j_max));
            T_total = Scalar{2} * T_j;
            dist = sign_v * abs_v * T_j;
            a_reached = false;
        }
        return {dist, T_j, T_total, a_reached};
    }

    Scalar compute_decel_distance(Scalar abs_v, Scalar T_j, Scalar T_const) const
    {
        auto const v1 = abs_v - this->m_j_max * T_j * T_j / Scalar{2};
        auto const q1 = abs_v * T_j - this->m_j_max * T_j * T_j * T_j / Scalar{6};
        auto const v2 = v1 - this->m_a_max * T_const;
        auto const q2 = q1 + v1 * T_const - this->m_a_max * T_const * T_const / Scalar{2};
        auto const q3 = q2 + v2 * T_j - this->m_a_max * T_j * T_j / Scalar{2} + this->m_j_max * T_j * T_j * T_j / Scalar{6};

        return q3;
    }

    void append_brake_phases(Scalar v0, stop_result const& info)
    {
        auto const sign_v = (v0 >= Scalar{0}) ? Scalar{1} : Scalar{-1};
        auto const first = this->m_phase_count;
        this->append_phase(info.T_j, -sign_v * this->m_j_max);
        if(info.a_max_reached)
        {
            this->append_phase(info.brake_duration - Scalar{2} * info.T_j, Scalar{0});
        }
        this->append_phase(info.T_j, sign_v * this->m_j_max);
        this->pin_group_end(first, Scalar{0}, Scalar{0});
    }
};

}

#endif
