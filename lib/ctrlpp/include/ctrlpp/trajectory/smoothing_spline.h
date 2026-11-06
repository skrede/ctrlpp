#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_SMOOTHING_SPLINE_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_SMOOTHING_SPLINE_H

// Cubic smoothing spline: (R + lambda Q^T Q) d = Q^T q, with lambda = 2 (1 - mu) / (3 mu), is solved for the second
// derivatives at the interior knots, d_0 = d_n = 0, and the coefficients follow from the smoothed positions q - lambda Q d.
// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
// Automatic Machines and Robots", 2008, Sec. 4.4.5
// @cite reinsch1967 -- Reinsch, "Smoothing by Spline Functions", Numerische Mathematik 10:177-183, 1967 (original smoothing-spline derivation)
// @cite deboor2001 -- de Boor, "A Practical Guide to Splines", Springer, 2001 (cubic smoothing splines)

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/trajectory_segment.h"

#include "ctrlpp/trajectory/detail/spline_checks.h"
#include "ctrlpp/trajectory/detail/smoothing_spline_solve.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <vector>
#include <cstddef>
#include <algorithm>

namespace ctrlpp
{

template <ctrlpp_floating_scalar Scalar>
class smoothing_spline
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
        Scalar mu{Scalar{0.5}};        // in (0, 1]: 1 interpolates, near 0 is smoothest
    };

    // The system's entries are bounded before the solve, by the square root of the largest Scalar, since the decomposition
    // squares them; on the smallest system an overflowed entry gives a finite, plausible and wrong answer rather than a NaN.
    // @cite biagiotti2009 -- Sec. 4.4.5
    static ctrlpp::expected<smoothing_spline, spline_error> create(config const& cfg)
    {
        if(auto const knots = detail::check_spline_knots(cfg.times, cfg.positions, Scalar{0}, Scalar{0}); !knots)
        {
            return ctrlpp::unexpected(knots.error());
        }
        // Written as a negated conjunction so a NaN mu fails both comparisons and is rejected.
        if(!(cfg.mu > Scalar{0} && cfg.mu <= Scalar{1}))
        {
            return ctrlpp::unexpected(spline_error::mu_out_of_range);
        }
        if(!detail::regularized_system_representable(cfg.times, cfg.mu))
        {
            return ctrlpp::unexpected(spline_error::unrepresentable_spline);
        }

        smoothing_spline spline{unchecked_t{}, cfg};
        if(!spline.representable(cfg))
        {
            return ctrlpp::unexpected(spline_error::unrepresentable_spline);
        }
        return spline;
    }

    // @cite biagiotti2009 -- Sec. 4.4.5
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
    bool m_coefficients_representable;

    // @cite biagiotti2009 -- Sec. 4.4.5
    smoothing_spline(unchecked_t, config const& cfg)
        : m_times{cfg.times}
        , m_coefficients_representable{true}
    {
        auto const n_pts = m_times.size();
        auto const n = n_pts - 1;

        auto const lambda = detail::smoothing_weight(cfg.mu);

        std::vector<Scalar> h(n);
        for(std::size_t i = 0; i < n; ++i)
        {
            h[i] = m_times[i + 1] - m_times[i];
        }

        // Two waypoints leave no interior knot: the spline is the line through them.
        if(n_pts == 2)
        {
            m_coefficients_representable = detail::build_linear_coefficients(m_coeffs, cfg.positions, h);
            return;
        }

        auto const m = n - 1;
        solve_and_build_coeffs(cfg.positions, h, lambda, n_pts, n, m);
    }

    void solve_and_build_coeffs(std::vector<Scalar> const& pos, std::vector<Scalar> const& h, Scalar lambda, std::size_t n_pts, std::size_t n, std::size_t m)
    {
        auto const mi = static_cast<Eigen::Index>(m);
        auto const ni = static_cast<Eigen::Index>(n_pts);

        detail::smoothing_matrix<Scalar> const R = detail::smoothness_matrix(h, mi);
        detail::smoothing_matrix<Scalar> const Q = detail::second_difference_matrix(h, ni, mi);
        detail::smoothing_matrix<Scalar> const QtQ = Q.transpose() * Q;
        detail::smoothing_matrix<Scalar> const A = R + lambda * QtQ;
        detail::smoothing_vector<Scalar> const q = detail::knot_vector(pos);
        detail::smoothing_vector<Scalar> const rhs = Q.transpose() * q;
        detail::smoothing_vector<Scalar> const d_int = A.colPivHouseholderQr().solve(rhs);

        std::vector<Scalar> d(n_pts, Scalar{0});
        for(Eigen::Index i = 0; i < mi; ++i)
        {
            d[static_cast<std::size_t>(i) + 1] = d_int(i);
        }
        detail::smoothing_vector<Scalar> const s_vec = q - lambda * Q * d_int;
        std::vector<Scalar> s(n_pts);
        for(Eigen::Index i = 0; i < ni; ++i)
        {
            s[static_cast<std::size_t>(i)] = s_vec(i);
        }
        m_coefficients_representable = detail::build_smoothing_coefficients(m_coeffs, pos, h, d, s, n);
    }

    bool representable(config const& cfg) const
    {
        if(!m_coefficients_representable || !std::isfinite(duration()))
        {
            return false;
        }
        auto const start = evaluate(cfg.times.front());
        auto const end = evaluate(cfg.times.back());
        return start.position.allFinite() && start.velocity.allFinite() && start.acceleration.allFinite() && end.position.allFinite() && end.velocity.allFinite() &&
                end.acceleration.allFinite();
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

static_assert(trajectory_segment<smoothing_spline<double>, double, 1>);

}

#endif
