#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_TRAPEZOIDAL_RESCALE_SHAPES_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_TRAPEZOIDAL_RESCALE_SHAPES_H

#include "ctrlpp/trajectory/detail/trapezoidal_state.h"

#include "ctrlpp/trajectory/trajectory_error.h"

#include "ctrlpp/expected.h"

#include <cmath>
#include <limits>
#include <cstdint>
#include <algorithm>

namespace ctrlpp::detail
{

enum class rescale_shape
{
    plateau,
    ramp_through,
    valley
};

template <typename Scalar>
struct duration_with_scale
{
    duration_with_scale()
        : T{}
        , scale{}
    {
    }

    Scalar T;
    Scalar scale;
};

// The solved shapes carry the distance from their own boundary as the unknown, never the cruise velocity
// itself: that distance is routinely orders of magnitude smaller than the boundary velocity it sits on.
template <typename Scalar>
class trapezoidal_rescale_shapes : protected trapezoidal_state<Scalar>
{
protected:
    using trapezoidal_state<Scalar>::trapezoidal_state;

    // Signed: negative when the smaller boundary velocity is the larger in magnitude, the case of an axis
    // moving away from its target faster than toward it, which the |v0^2 - v1^2| spelling loses.
    Scalar ramp_through_distance() const
    {
        auto const v_lo = std::min(this->m_v0, this->m_v1);
        auto const v_hi = std::max(this->m_v0, this->m_v1);
        return (v_hi - v_lo) * (v_hi + v_lo) / (Scalar{2} * this->m_a);
    }

    // The scale is the largest operand that entered, in units of time, not the duration alone: the cruise
    // numerator cancels to nothing whenever the two ramps sweep almost all of the displacement.
    duration_with_scale<Scalar> duration_and_scale_at(Scalar v) const
    {
        auto const h = std::abs(this->m_q1 - this->m_q0);
        auto const T_a = std::abs(v - this->m_v0) / this->m_a;
        auto const T_d = std::abs(v - this->m_v1) / this->m_a;
        auto const d_a = (this->m_v0 + v) * T_a / Scalar{2};
        auto const d_d = (this->m_v1 + v) * T_d / Scalar{2};
        auto const T = T_a + (h - d_a - d_d) / v + T_d;
        auto const cruise_scale = std::max({h, std::abs(d_a), std::abs(d_d)}) / std::abs(v);
        duration_with_scale<Scalar> result;
        result.T = T;
        result.scale = std::max({T_a, T_d, cruise_scale, std::abs(T)});
        return result;
    }

    // The ramp-through residual collapses whenever the acceleration was raised, since the raise makes the ramp
    // distance equal the displacement. Six chained roundings form it, each worth one unit in the last place at
    // the larger of h and (v0^2 + v1^2) / (2a).
    static bool residual_resolved(Scalar ramp_residual, Scalar residual_scale)
    {
        constexpr std::int32_t residual_rounding_ops = 6;
        auto const residual_floor = Scalar{residual_rounding_ops} * std::numeric_limits<Scalar>::epsilon() * residual_scale;
        return ramp_residual > residual_floor;
    }

    // duration_and_scale_at chains fifteen arithmetic operations, each counted whether or not it rounds, so the
    // count bounds the boundary duration's error from above. The floors are tested before any sign gate,
    // because a coefficient built from an unresolved residual carries a sign that is noise.
    static bool boundary_resolved(bool ramp_resolved, Scalar dT, Scalar boundary_scale)
    {
        constexpr std::int32_t boundary_duration_rounding_ops = 15;
        auto const boundary_floor = Scalar{boundary_duration_rounding_ops} * std::numeric_limits<Scalar>::epsilon() * boundary_scale;
        return ramp_resolved && dT > boundary_floor;
    }

    bool peak_below_boundary(Scalar h, Scalar a, Scalar v_sum_sq, Scalar v_hi) const
    {
        auto const v_tri = std::sqrt(a * h + v_sum_sq / Scalar{2});
        auto const v_top = std::min(v_tri, this->m_v_max);
        return v_top < v_hi;
    }

    // The shape that reaches a vanishing cruise velocity fixes the supremum of the reachable durations, which is
    // infinite exactly when that shape's residual displacement is positive. With both boundary velocities
    // positive the valley decides it at its own solve, where the test is exact.
    bool beyond_supremum(Scalar h, Scalar a, Scalar v_lo, Scalar v_hi, bool ramp_resolved, Scalar T_new) const
    {
        bool unbounded{true};
        Scalar T_sup{};
        if(v_hi <= Scalar{0})
        {
            unbounded = (h > Scalar{0});
            T_sup = -(this->m_v0 + this->m_v1) / a;
        }
        else if(v_lo <= Scalar{0})
        {
            unbounded = ramp_resolved;
            T_sup = std::abs(this->m_v1 - this->m_v0) / a;
        }
        return !unbounded && T_new > T_sup;
    }

    // v = v_hi + e gives e^2 - B e + C = 0 with B = a (r / v_hi - dT) and C = a v_hi dT; the smaller positive
    // root, in the form that adds. B is a difference of two terms that agree to the full significand deep in
    // the plateau: six roundings form r, one divides it, fifteen form the boundary duration, one subtracts T_new.
    ctrlpp::expected<Scalar, trajectory_error> plateau_rise_root(Scalar a, Scalar v_hi, Scalar residual_scale, Scalar ramp_residual, Scalar boundary_scale, Scalar dT) const
    {
        constexpr std::int32_t coefficient_rounding_ops = 23;
        auto const coefficient_scale = std::max(residual_scale / v_hi, boundary_scale);
        auto const coefficient_floor = Scalar{coefficient_rounding_ops} * std::numeric_limits<Scalar>::epsilon() * coefficient_scale;
        auto const rise_coefficient = ramp_residual / v_hi - dT;
        if(!(std::abs(rise_coefficient) > coefficient_floor))
        {
            return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
        }
        auto const B = a * rise_coefficient;
        auto const C = a * v_hi * dT;
        if(!(B > Scalar{0}))
        {
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }
        auto const disc = B * B - Scalar{4} * C;
        if(!std::isfinite(disc))
        {
            return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
        }
        if(!(disc >= Scalar{0}))
        {
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }
        return v_hi + Scalar{2} * C / (B + std::sqrt(disc));
    }

    ctrlpp::expected<Scalar, trajectory_error> plateau_rise_velocity(
            Scalar a, Scalar v_hi, Scalar residual_scale, Scalar ramp_residual, bool ramp_resolved, duration_with_scale<Scalar> const& at_hi, Scalar T_new) const
    {
        auto const dT = at_hi.T - T_new;
        if(!boundary_resolved(ramp_resolved, dT, at_hi.scale))
        {
            return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
        }
        return plateau_rise_root(a, v_hi, residual_scale, ramp_residual, at_hi.scale, dT);
    }

    ctrlpp::expected<Scalar, trajectory_error> ramp_through_velocity(Scalar a, Scalar ramp_residual, bool ramp_resolved, Scalar T_new) const
    {
        auto const denom = T_new - std::abs(this->m_v1 - this->m_v0) / a;
        if(!ramp_resolved || !(denom > Scalar{0}))
        {
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }
        return ramp_residual / denom;
    }
};

}

#endif
