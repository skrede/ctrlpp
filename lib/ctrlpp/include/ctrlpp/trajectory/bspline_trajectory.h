#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_BSPLINE_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_BSPLINE_TRAJECTORY_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
// Automatic Machines and Robots", 2008, Sec. 4.5
// @cite deboor2001 -- de Boor, "A Practical Guide to Splines", Springer, 2001 (de Boor evaluation algorithm)
// @cite piegl1997 -- Piegl & Tiller, "The NURBS Book", 2nd ed., 1997, Ch. 2-3 (B-spline basis and knot averaging)

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/trajectory_segment.h"

#include "ctrlpp/trajectory/detail/bspline_basis.h"
#include "ctrlpp/trajectory/detail/bspline_checks.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include <Eigen/Dense>

#include <vector>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <algorithm>

namespace ctrlpp
{

// @cite biagiotti2009 -- Sec. 4.5
template <ctrlpp_floating_scalar Scalar, int Degree>
class bspline_trajectory
{
    static_assert(Degree >= 1, "bspline_trajectory requires Degree >= 1");

    struct unchecked_t
    {
        explicit unchecked_t() = default;
    };

public:
    struct config
    {
        std::vector<Scalar> control_points; // n+1 control points
        std::vector<Scalar> knot_vector{};  // empty: a uniform clamped knot vector is generated
    };

    static ctrlpp::expected<bspline_trajectory, spline_error> create(config const& cfg)
    {
        if(cfg.control_points.size() < static_cast<std::size_t>(Degree) + 1)
            return ctrlpp::unexpected(spline_error::too_few_control_points);
        if(!detail::all_finite(cfg.control_points) || !detail::all_finite(cfg.knot_vector))
            return ctrlpp::unexpected(spline_error::non_finite_input);
        if(auto const knots = detail::check_knots<Degree>(cfg.control_points.size(), cfg.knot_vector); !knots)
            return ctrlpp::unexpected(knots.error());

        bspline_trajectory spline{unchecked_t{}, cfg};
        if(!detail::all_finite(spline.m_d1_points) || !detail::all_finite(spline.m_d2_points))
            return ctrlpp::unexpected(spline_error::unrepresentable_spline);
        return spline;
    }

    trajectory_point<Scalar, 1> evaluate(Scalar t) const
    {
        constexpr std::int32_t p = Degree;
        auto const n = static_cast<std::int32_t>(m_control_points.size()) - 1;

        auto const t_min = m_knots[static_cast<std::size_t>(p)];
        auto const t_max = m_knots[static_cast<std::size_t>(n + 1)];
        auto const tc = std::clamp(t, t_min, t_max);

        auto const pos = detail::de_boor<Degree>(tc, m_control_points, m_knots, p);
        auto const vel = detail::de_boor<Degree>(tc, m_d1_points, m_d1_knots, p - 1);
        auto acc = Scalar{0};
        if constexpr(p >= 2)
            acc = detail::de_boor<Degree>(tc, m_d2_points, m_d2_knots, p - 2);

        return {
                .position = Vector<Scalar, 1>{pos},
                .velocity = Vector<Scalar, 1>{vel},
                .acceleration = Vector<Scalar, 1>{acc},
        };
    }

    Scalar duration() const
    {
        constexpr std::int32_t p = Degree;
        auto const n = static_cast<std::int32_t>(m_control_points.size()) - 1;
        return m_knots[static_cast<std::size_t>(n + 1)] - m_knots[static_cast<std::size_t>(p)];
    }

private:
    std::vector<Scalar> m_control_points;
    std::vector<Scalar> m_knots;
    std::vector<Scalar> m_d1_points;
    std::vector<Scalar> m_d1_knots;
    std::vector<Scalar> m_d2_points;
    std::vector<Scalar> m_d2_knots;

    bspline_trajectory(unchecked_t, config const& cfg)
        : m_control_points{cfg.control_points}
    {
        auto const n = static_cast<std::int32_t>(m_control_points.size()) - 1;
        constexpr std::int32_t p = Degree;

        if(cfg.knot_vector.empty())
        {
            auto const m = n + p + 1;
            m_knots.resize(static_cast<std::size_t>(m + 1));
            fill_uniform_knots(n, m);
        }
        else
            m_knots = cfg.knot_vector;

        // The derivative data is formed once here so evaluate stays allocation free.
        auto d1 = derivative_data(m_control_points, m_knots, p);
        m_d1_points = std::move(d1.first);
        m_d1_knots = std::move(d1.second);

        if constexpr(p >= 2)
        {
            auto d2 = derivative_data(m_d1_points, m_d1_knots, p - 1);
            m_d2_points = std::move(d2.first);
            m_d2_knots = std::move(d2.second);
        }
    }

    // The uniform clamped knot vector: p + 1 zeros, uniform interior knots, p + 1 ones.
    void fill_uniform_knots(std::int32_t n, std::int32_t m)
    {
        constexpr std::int32_t p = Degree;
        for(std::int32_t i = 0; i <= p; ++i)
            m_knots[static_cast<std::size_t>(i)] = Scalar{0};
        for(std::int32_t i = m - p; i <= m; ++i)
            m_knots[static_cast<std::size_t>(i)] = Scalar{1};
        auto const num_interior = n - p;
        for(std::int32_t i = 1; i <= num_interior; ++i)
            m_knots[static_cast<std::size_t>(p + i)] = static_cast<Scalar>(i) / static_cast<Scalar>(num_interior + 1);
    }

    // Derivative control points Q_i = deg (P_{i+1} - P_i) / (U_{i+deg+1} - U_{i+1}) over the knot vector without its end knots.
    // @cite biagiotti2009 -- Sec. 4.5; App. B.1.4, eq. (B.8)-(B.10), p.476
    static std::pair<std::vector<Scalar>, std::vector<Scalar>> derivative_data(std::vector<Scalar> const& points, std::vector<Scalar> const& U, std::int32_t deg)
    {
        auto const n = static_cast<std::int32_t>(points.size()) - 1;

        std::vector<Scalar> d_points(static_cast<std::size_t>(n));
        for(std::int32_t i = 0; i < n; ++i)
        {
            auto const denom = U[static_cast<std::size_t>(i + deg + 1)] - U[static_cast<std::size_t>(i + 1)];
            auto const rise = points[static_cast<std::size_t>(i + 1)] - points[static_cast<std::size_t>(i)];
            // 0/0 := 0, Piegl & Tiller Eq. 2.5.
            d_points[static_cast<std::size_t>(i)] = denom == Scalar{0} ? Scalar{0} : static_cast<Scalar>(deg) * rise / denom;
        }

        std::vector<Scalar> d_knots(U.begin() + 1, U.end() - 1);

        return {d_points, d_knots};
    }
};

template <typename Scalar>
Scalar basis_function(std::int32_t i, std::int32_t p, Scalar t, std::vector<Scalar> const& U)
{
    return detail::bspline_basis(i, p, t, U);
}

template <typename Scalar, int Degree>
ctrlpp::expected<bspline_trajectory<Scalar, Degree>, spline_error> make_bspline_interpolation(std::vector<Scalar> const& times, std::vector<Scalar> const& positions)
{
    if(times.size() != positions.size())
        return ctrlpp::unexpected(spline_error::size_mismatch);
    if(times.size() < static_cast<std::size_t>(Degree + 1))
        return ctrlpp::unexpected(spline_error::too_few_points);
    if(auto const valid = detail::check_interpolation_input(times, positions); !valid)
        return ctrlpp::unexpected(valid.error());

    auto knots = detail::averaged_knots<Degree>(times);
    auto const N = detail::collocation_matrix<Degree>(times, knots);

    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> const rhs = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1> const>(positions.data(), N.rows());
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> ctrl = N.colPivHouseholderQr().solve(rhs);

    std::vector<Scalar> control_points(ctrl.data(), ctrl.data() + ctrl.size());
    if(!detail::all_finite(control_points))
        return ctrlpp::unexpected(spline_error::unrepresentable_spline);

    return bspline_trajectory<Scalar, Degree>::create({
            .control_points = std::move(control_points),
            .knot_vector = std::move(knots),
    });
}

static_assert(trajectory_segment<bspline_trajectory<double, 3>, double, 1>);
static_assert(trajectory_segment<bspline_trajectory<double, 5>, double, 1>);

}

#endif
