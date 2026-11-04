#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_DOUBLE_S_SOLVE_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_DOUBLE_S_SOLVE_H

#include "ctrlpp/trajectory/detail/jerk_limited_ramp.h"
#include "ctrlpp/trajectory/detail/double_s_rest_to_rest.h"

#include <cmath>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar>
class double_s_solve : protected double_s_rest_to_rest<Scalar>
{
    using base = double_s_rest_to_rest<Scalar>;

protected:
    // The zero-displacement branch is keyed on exact equality: a displacement small enough to underflow the
    // rest of the algebra is a rejection `create` makes, not a command rounded down to a standstill.
    // @cite biagiotti2009 -- Sec. 3.4, eq. (3.17)-(3.33), Fig. 3.18, p.79-96
    template <typename Config>
    explicit double_s_solve(Config const& cfg)
        : base{cfg}
    {
        auto const h_signed = cfg.q1 - cfg.q0;
        if(h_signed == Scalar{0})
        {
            this->m_sigma = Scalar{1};
            this->m_degenerate = true;
            return;
        }

        // @cite biagiotti2009 -- Sec. 3.4.2, eq. (3.31)-(3.33), p.87
        this->m_sigma = (h_signed > Scalar{0}) ? Scalar{1} : Scalar{-1};
        auto const h = this->m_sigma * h_signed;
        auto const v_max = cfg.v_max;
        auto const a_max = cfg.a_max;
        this->m_pv0 = this->m_sigma * cfg.v0;
        this->m_pv1 = this->m_sigma * cfg.v1;
        if(cfg.v0 == Scalar{0} && cfg.v1 == Scalar{0})
        {
            this->compute_zero_bc(h, v_max, a_max, cfg.j_max);
        }
        else
        {
            compute_general_bc(h, v_max, a_max, cfg.j_max);
        }
    }

    // In the positive-displacement frame; reduces exactly to the zero-boundary case when pv0 = pv1 = 0.
    // @cite biagiotti2009 -- Sec. 3.4.1, eq. (3.19)-(3.27), Fig. 3.18, p.79-85
    void compute_general_bc(Scalar h, Scalar v_max, Scalar a_max, Scalar j_max)
    {
        auto const pv_hi = std::max(this->m_pv0, this->m_pv1);
        auto const pv_lo = std::min(this->m_pv0, this->m_pv1);

        // A cruise at v_max exists only when the limit is at least both boundary velocities; otherwise the
        // peak lies above both of them and only the cruise-free construction can express that.
        if(v_max >= pv_hi && set_general_cruise(h, v_max, a_max, j_max))
        {
            return;
        }

        auto const r_max = v_max - pv_hi;
        auto const rise = general_rise(h, pv_hi, pv_lo, r_max, a_max, j_max);
        // Formed from the nearer of v_max and pv_hi, the peak is at most v_max: r_max - rise >= 0 on one side and
        // pv_hi + rise <= (v_max + pv_hi) / 2 to within rounding on the other.
        auto const v_peak = (rise + rise > r_max) ? v_max - (r_max - rise) : pv_hi + rise;
        set_general_no_cruise(rise, pv_hi, pv_lo, v_peak, a_max, j_max);
    }

    // @cite biagiotti2009 -- Sec. 3.4.1, eq. (3.19)-(3.21), p.80-81
    bool set_general_cruise(Scalar h, Scalar v_max, Scalar a_max, Scalar j_max)
    {
        auto const accel = make_velocity_ramp(v_max - this->m_pv0, a_max, j_max);
        auto const decel = make_velocity_ramp(v_max - this->m_pv1, a_max, j_max);
        auto const swept = (this->m_pv0 + v_max) * accel.T / Scalar{2} + (v_max + this->m_pv1) * decel.T / Scalar{2};
        auto const T_v = (h - swept) / v_max;

        if(!(T_v > Scalar{0}))
        {
            return false;
        }
        this->set_ramps(accel, decel);
        this->m_T_v = T_v;
        this->m_v_lim = v_max;
        this->m_degenerate = (accel.a_lim < a_max) || (decel.a_lim < a_max);
        this->m_T = this->m_T_a + this->m_T_v + this->m_T_d;
        return true;
    }

    // The unknown is the peak's rise above pv_hi, not the peak: the rise is routinely orders of magnitude
    // below pv_hi, and subtracting pv_hi back out of a stored peak would discard most of its digits.
    // Where both ramps reach a_max the rise x solves x^2 + linear x + constant = 0; the larger root is
    // taken in the form that avoids subtracting nearly equal terms.
    // @cite biagiotti2009 -- Sec. 3.4.1, p.83-84
    static Scalar general_rise(Scalar h, Scalar pv_hi, Scalar pv_lo, Scalar r_max, Scalar a_max, Scalar j_max)
    {
        auto const a_sq_over_j = a_max * a_max / j_max;
        auto const floor_rise = std::max(Scalar{0}, -pv_hi);

        Scalar rise{};
        if(h >= no_cruise_displacement(a_sq_over_j, pv_hi, pv_lo, a_max, j_max))
        {
            auto const linear = Scalar{2} * pv_hi + a_sq_over_j;
            auto const constant = (pv_hi - pv_lo) * (pv_hi + pv_lo) / Scalar{2} + a_sq_over_j * (Scalar{3} * pv_hi + pv_lo) / Scalar{2} - a_max * h;
            auto const root = std::sqrt(linear * linear - Scalar{4} * constant);
            rise = (linear >= Scalar{0}) ? (Scalar{-2} * constant / (linear + root)) : ((root - linear) / Scalar{2});
        }
        else
        {
            rise = solve_no_cruise_rise(h, pv_hi, pv_lo, a_max, j_max, floor_rise, a_sq_over_j);
        }
        // The cruise-free branch has no rise above r_max, so a rise beyond it is re-solved on that domain.
        if(rise > r_max)
        {
            rise = solve_no_cruise_rise(h, pv_hi, pv_lo, a_max, j_max, floor_rise, r_max);
        }
        return rise;
    }

    void set_general_no_cruise(Scalar rise, Scalar pv_hi, Scalar pv_lo, Scalar v_peak, Scalar a_max, Scalar j_max)
    {
        auto const near_ramp = make_velocity_ramp(rise, a_max, j_max);
        auto const far_ramp = make_velocity_ramp(rise + (pv_hi - pv_lo), a_max, j_max);
        auto const& accel = (this->m_pv0 >= this->m_pv1) ? near_ramp : far_ramp;
        auto const& decel = (this->m_pv0 >= this->m_pv1) ? far_ramp : near_ramp;

        this->set_ramps(accel, decel);
        this->m_T_v = Scalar{0};
        this->m_v_lim = v_peak;
        this->m_degenerate = true;
        this->m_T = this->m_T_a + this->m_T_v + this->m_T_d;
    }

    // Bisection on the rise of the cruise-free profile that sweeps exactly h, where the two ramps carry
    // different square roots of the rise and no elementary closed form exists. The bracket starts at a zero
    // peak, above which the swept distance is monotone in the rise; it ends when the rise bracket is exhausted
    // or the swept-distance bracket stops shrinking.
    static Scalar solve_no_cruise_rise(Scalar h, Scalar pv_hi, Scalar pv_lo, Scalar a_max, Scalar j_max, Scalar lo, Scalar hi)
    {
        auto swept_lo = no_cruise_displacement(lo, pv_hi, pv_lo, a_max, j_max);
        auto swept_hi = no_cruise_displacement(hi, pv_hi, pv_lo, a_max, j_max);
        auto width = swept_hi - swept_lo;

        while(halve_rise_bracket(h, pv_hi, pv_lo, a_max, j_max, lo, hi, swept_lo, swept_hi))
        {
            auto const next_width = swept_hi - swept_lo;
            if(!(next_width < width))
            {
                break;
            }
            width = next_width;
        }

        return (std::abs(swept_lo - h) <= std::abs(swept_hi - h)) ? lo : hi;
    }

    static bool halve_rise_bracket(Scalar h, Scalar pv_hi, Scalar pv_lo, Scalar a_max, Scalar j_max, Scalar& lo, Scalar& hi, Scalar& swept_lo, Scalar& swept_hi)
    {
        auto const mid = lo + (hi - lo) / Scalar{2};
        if(!(mid > lo) || !(mid < hi))
        {
            return false;
        }
        auto const swept_mid = no_cruise_displacement(mid, pv_hi, pv_lo, a_max, j_max);
        if(swept_mid < h)
        {
            lo = mid;
            swept_lo = swept_mid;
        }
        else
        {
            hi = mid;
            swept_hi = swept_mid;
        }
        return true;
    }
};

}

#endif
