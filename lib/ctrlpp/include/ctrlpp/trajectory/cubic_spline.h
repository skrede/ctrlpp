#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_CUBIC_SPLINE_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_CUBIC_SPLINE_H

// C2 cubic spline interpolation in the velocity formulation: one tridiagonal system for the knot velocities.
// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
// Automatic Machines and Robots", 2009, Sec. 4.4, eq. (4.10)-(4.11)
// @cite deboor2001 -- de Boor, "A Practical Guide to Splines", Springer, 2001 (canonical reference for cubic-spline interpolation)

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/trajectory_segment.h"

#include "ctrlpp/trajectory/detail/spline_checks.h"
#include "ctrlpp/trajectory/detail/cubic_spline_solve.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include <array>
#include <cmath>
#include <limits>
#include <vector>
#include <cstddef>
#include <algorithm>

namespace ctrlpp
{

enum class boundary_condition
{
    natural, // zero acceleration at the ends (M_0 = M_n = 0)
    clamped, // assigned end velocities v0 and vn
    periodic // v_0 = v_n and a_0 = a_n, requiring q_0 = q_n
};

template <ctrlpp_floating_scalar Scalar>
class cubic_spline
{
    struct unchecked_t
    {
        explicit unchecked_t() = default;
    };

public:
    struct config
    {
        std::vector<Scalar> times;     // t_0 ... t_n
        std::vector<Scalar> positions; // q_0 ... q_n
        boundary_condition bc{boundary_condition::natural};
        Scalar v0{}; // clamped only
        Scalar vn{}; // clamped only
    };

    // The realized coefficients are checked rather than predicted from the knot spacing: the cubic coefficient runs with the
    // waypoint separation over the cube of the span, and its numerator carries a waypoint scale that collinear data cancels.
    // @cite biagiotti2009 -- Sec. 4.4, eq. (4.10)-(4.11)
    static ctrlpp::expected<cubic_spline, spline_error> create(config const& cfg)
    {
        if(auto const checked = check_config(cfg); !checked)
        {
            return ctrlpp::unexpected(checked.error());
        }

        cubic_spline spline{unchecked_t{}, cfg};
        if(!spline.m_system_solved || !spline.m_coefficients_representable)
        {
            return ctrlpp::unexpected(spline_error::unrepresentable_spline);
        }
        return spline;
    }

    // @cite biagiotti2009 -- Sec. 4.4
    trajectory_point<Scalar, 1> evaluate(Scalar t) const
    {
        auto const tc = std::clamp(t, m_times.front(), m_times.back());
        auto const span = find_span(tc);

        auto const s = tc - m_times[span];
        auto const& c = m_coeffs[span];

        auto const q = c[0] + s * (c[1] + s * (c[2] + s * c[3]));
        auto const dq = c[1] + s * (Scalar{2} * c[2] + Scalar{3} * c[3] * s);
        auto const ddq = Scalar{2} * c[2] + Scalar{6} * c[3] * s;

        return {
                .position = Vector<Scalar, 1>{q},
                .velocity = Vector<Scalar, 1>{dq},
                .acceleration = Vector<Scalar, 1>{ddq},
        };
    }

    Scalar duration() const { return m_times.back() - m_times.front(); }

private:
    std::vector<Scalar> m_times;
    std::vector<std::array<Scalar, 4>> m_coeffs; // {a, b, c, d} of a + b s + c s^2 + d s^3, s = t - t_i
    bool m_system_solved;
    bool m_coefficients_representable;

    // @cite biagiotti2009 -- Sec. 4.4, eq. (4.10)-(4.11)
    cubic_spline(unchecked_t, config const& cfg)
        : m_times{cfg.times}
        , m_system_solved{true}
        , m_coefficients_representable{true}
    {
        auto const n_pts = m_times.size();
        auto const n = n_pts - 1;
        std::vector<Scalar> h(n);
        std::vector<Scalar> delta(n);
        for(std::size_t i = 0; i < n; ++i)
        {
            h[i] = m_times[i + 1] - m_times[i];
            delta[i] = (cfg.positions[i + 1] - cfg.positions[i]) / h[i];
        }
        std::vector<Scalar> vel(n_pts);
        auto const solved = (cfg.bc == boundary_condition::natural) ? detail::solve_natural(h, delta, cfg.positions, vel, n_pts, n)
                : (cfg.bc == boundary_condition::clamped)           ? detail::solve_clamped(h, delta, cfg.positions, vel, n_pts, n, cfg.v0, cfg.vn)
                                                                    : detail::solve_periodic(h, delta, cfg.positions, vel, n_pts, n);
        if(!solved)
        {
            m_system_solved = false;
            return;
        }
        auto const position_scale = detail::largest_magnitude(cfg.positions);
        m_coeffs.resize(n);
        for(std::size_t i = 0; i < n; ++i)
        {
            set_span(i, cfg.positions, h, vel, position_scale);
        }
    }

    void set_span(std::size_t i, std::vector<Scalar> const& pos, std::vector<Scalar> const& h, std::vector<Scalar> const& vel, Scalar position_scale)
    {
        auto const qi = pos[i];
        auto const qi1 = pos[i + 1];
        auto const vi = vel[i];
        auto const vi1 = vel[i + 1];
        auto const hi = h[i];

        auto const quadratic_numerator = Scalar{3} * (qi1 - qi) / hi - Scalar{2} * vi - vi1;
        auto const cubic_numerator = Scalar{2} * (qi - qi1) / hi + vi + vi1;

        m_coeffs[i][0] = qi;
        m_coeffs[i][1] = vi;
        m_coeffs[i][2] = quadratic_numerator / hi;
        m_coeffs[i][3] = cubic_numerator / (hi * hi);

        m_coefficients_representable = m_coefficients_representable && detail::cubic_coefficients_representable(m_coeffs[i], quadratic_numerator, cubic_numerator, hi, position_scale);
    }

    static ctrlpp::expected<void, spline_error> check_config(config const& cfg)
    {
        if(auto const knots = detail::check_spline_knots(cfg.times, cfg.positions, cfg.v0, cfg.vn); !knots || cfg.bc != boundary_condition::periodic)
        {
            return knots;
        }
        if(cfg.times.size() < 3)
        {
            return ctrlpp::unexpected(spline_error::periodic_too_few_points);
        }
        // Relative to the endpoint magnitude, so it holds across scales and precisions; the margin is the difference's rounding budget.
        constexpr Scalar periodic_match_ulps = Scalar{4};
        auto const periodic_scale = std::max({Scalar{1}, std::abs(cfg.positions.front()), std::abs(cfg.positions.back())});
        auto const endpoint_difference = std::abs(cfg.positions.front() - cfg.positions.back());
        if(!std::isfinite(endpoint_difference) || endpoint_difference > periodic_match_ulps * std::numeric_limits<Scalar>::epsilon() * periodic_scale)
        {
            return ctrlpp::unexpected(spline_error::periodic_endpoint_mismatch);
        }
        return {};
    }

    std::size_t find_span(Scalar t) const
    {
        auto const n_spans = m_coeffs.size();
        auto it = std::upper_bound(m_times.begin(), m_times.end(), t);
        if(it == m_times.begin())
        {
            return 0;
        }
        auto idx = static_cast<std::size_t>(it - m_times.begin()) - 1;
        if(idx >= n_spans)
        {
            idx = n_spans - 1;
        }
        return idx;
    }
};

static_assert(trajectory_segment<cubic_spline<double>, double, 1>);

}

#endif
