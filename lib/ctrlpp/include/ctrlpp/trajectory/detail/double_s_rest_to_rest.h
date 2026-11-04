#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_DOUBLE_S_REST_TO_REST_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_DOUBLE_S_REST_TO_REST_H

#include "ctrlpp/trajectory/detail/double_s_state.h"
#include "ctrlpp/trajectory/detail/jerk_limited_ramp.h"

namespace ctrlpp::detail
{

template <typename Scalar>
class double_s_rest_to_rest : protected double_s_state<Scalar>
{
protected:
    using double_s_state<Scalar>::double_s_state;

    // @cite biagiotti2009 -- Sec. 3.4.3, p.88-91
    void compute_zero_bc(Scalar h, Scalar v_max, Scalar a_max, Scalar j_max)
    {
        // @cite biagiotti2009 -- Sec. 3.4.3, eq. before (3.34), p.89
        bool const a_max_reached = (v_max * j_max >= a_max * a_max);
        auto const ramp = make_velocity_ramp(v_max, a_max, j_max);
        // The cruise is decided on the distance the cruise-free profile at v_max lacks, the c that forms its peak's deficit.
        auto const c = v_max * ramp.T - h;

        if(c < Scalar{0})
        {
            set_zero_bc_cruise(ramp, c, v_max, a_max_reached);
        }
        else
        {
            set_zero_bc_no_cruise(c, h, v_max, a_max, j_max, a_max_reached);
        }

        this->m_T = this->m_T_a + this->m_T_v + this->m_T_d;
    }

    void set_zero_bc_cruise(velocity_ramp<Scalar> const& ramp, Scalar c, Scalar v_max, bool a_max_reached)
    {
        this->m_T_j1 = ramp.T_j;
        this->m_T_j2 = ramp.T_j;
        this->m_a_lim_a = ramp.a_lim;
        this->m_a_lim_d = ramp.a_lim;
        this->m_T_a = ramp.T;
        this->m_T_d = ramp.T;
        this->m_T_v = -c / v_max;
        this->m_v_lim = v_max;
        this->m_degenerate = !a_max_reached;
    }

    // @cite biagiotti2009 -- Sec. 3.4.3, p.89-91
    void set_zero_bc_no_cruise(Scalar c, Scalar h, Scalar v_max, Scalar a_max, Scalar j_max, bool a_max_reached)
    {
        this->m_degenerate = true;
        this->m_T_v = Scalar{0};

        if(a_max_reached)
        {
            this->m_v_lim = no_cruise_peak(c, h, v_max, a_max, j_max);
            this->m_T_j1 = a_max / j_max;
            this->m_T_j2 = this->m_T_j1;
            this->m_T_a = this->m_T_j1 + this->m_v_lim / a_max;
            this->m_T_d = this->m_T_a;
            this->m_a_lim_a = a_max;
            this->m_a_lim_d = a_max;

            if(this->m_T_a < Scalar{2} * this->m_T_j1)
            {
                solve_doubly_degenerate(h, a_max, j_max);
            }
        }
        else
        {
            solve_doubly_degenerate(h, a_max, j_max);
        }
    }

    // @cite biagiotti2009 -- Sec. 3.4.3, p.90-91
    void solve_doubly_degenerate(Scalar h, Scalar a_max, Scalar j_max)
    {
        auto const T_j = doubly_degenerate_jerk_time(h, this->m_v_max, a_max, j_max);

        this->m_T_j1 = T_j;
        this->m_T_j2 = T_j;
        this->m_T_a = Scalar{2} * T_j;
        this->m_T_d = Scalar{2} * T_j;
        this->m_T_v = Scalar{0};
        this->m_a_lim_a = j_max * T_j;
        this->m_a_lim_d = this->m_a_lim_a;
        this->m_v_lim = this->m_a_lim_a * T_j;
    }
};

}

#endif
