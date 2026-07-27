#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_BSPLINE_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_BSPLINE_TRAJECTORY_H

/// @brief B-spline trajectory with compile-time degree and de Boor evaluation.
///
/// Supports configurable degree (cubic = 3, quintic = 5, etc.), auto-generated
/// uniform clamped knot vectors, and user-provided custom knot vectors. The
/// factory function make_bspline_interpolation() solves for control points that
/// pass through specified waypoints at given parameter values.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
/// Automatic Machines and Robots", 2009, Sec. 4.5
/// @cite deboor2001 -- de Boor, "A Practical Guide to Splines", Springer, 2001 (de Boor evaluation algorithm)
/// @cite piegl1997 -- Piegl & Tiller, "The NURBS Book", 2nd ed., 1997, Ch. 2-3 (B-spline basis and knot averaging)

#include "ctrlpp/expected.h"

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/trajectory_segment.h"

#include "ctrlpp/util/concepts.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <limits>
#include <vector>
#include <cstddef>
#include <utility>
#include <algorithm>

namespace ctrlpp
{

/// @brief B-spline trajectory with compile-time degree and de Boor evaluation.
///
/// Construction accepts control points and an optional knot vector. If no knot
/// vector is provided, a uniform clamped knot vector is generated automatically.
/// Evaluation uses de Boor's algorithm for numerically stable computation.
///
/// @tparam Scalar  Floating-point type (float, double, long double)
/// @tparam Degree  B-spline degree (compile-time, e.g. 3 for cubic, 5 for quintic)
///
/// @cite biagiotti2009 -- Sec. 4.5
template <ctrlpp_floating_scalar Scalar, int Degree>
class bspline_trajectory
{
    static_assert(Degree >= 1, "bspline_trajectory requires Degree >= 1");

  public:
    struct config
    {
        std::vector<Scalar> control_points; ///< n+1 control points
        std::vector<Scalar> knot_vector{};  ///< If empty, auto-generate uniform clamped
    };

    /// @brief Validate the configuration and construct a B-spline trajectory.
    ///
    /// Rejections, checked in order:
    ///  * fewer than Degree + 1 control points -> spline_error::too_few_control_points
    ///  * knot vector size differs from
    ///    control_points.size() + Degree + 1   -> spline_error::bad_knot_count
    ///  * knot vector not non-decreasing       -> spline_error::non_monotonic_knots
    ///
    /// An empty knot vector skips the knot checks; a uniform clamped knot vector
    /// is generated instead, which is valid by construction.
    ///
    /// @cite biagiotti2009 -- Sec. 4.5
    static auto create(config const& cfg)
        -> ctrlpp::expected<bspline_trajectory, spline_error>
    {
        auto const n = static_cast<int>(cfg.control_points.size()) - 1;
        constexpr int p = Degree;

        if (n < p) {
            return ctrlpp::unexpected(spline_error::too_few_control_points);
        }

        if (!cfg.knot_vector.empty()) {
            // Knot vector size must be n + p + 2
            auto const expected_size = static_cast<std::size_t>(n + p + 2);
            if (cfg.knot_vector.size() != expected_size) {
                return ctrlpp::unexpected(spline_error::bad_knot_count);
            }

            // Knots must be non-decreasing
            for (std::size_t i = 1; i < cfg.knot_vector.size(); ++i) {
                if (cfg.knot_vector[i] < cfg.knot_vector[i - 1]) {
                    return ctrlpp::unexpected(spline_error::non_monotonic_knots);
                }
            }
        }

        return bspline_trajectory{unchecked_t{}, cfg};
    }

    /// @brief Evaluate B-spline at parameter t, returning position, velocity, acceleration.
    ///
    /// Clamps t to the active parameter range [U[p], U[n+1]].
    ///
    /// @cite biagiotti2009 -- Sec. 4.5
    auto evaluate(Scalar t) const -> trajectory_point<Scalar, 1>
    {
        constexpr int p = Degree;
        auto const n = static_cast<int>(control_points_.size()) - 1;

        // Active parameter range
        auto const t_min = knots_[static_cast<std::size_t>(p)];
        auto const t_max = knots_[static_cast<std::size_t>(n + 1)];
        auto const tc = std::clamp(t, t_min, t_max);

        auto const pos = de_boor(tc, control_points_, knots_, p);

        // Velocity: derivative control points of degree p-1, precomputed at
        // construction since they depend only on the stored control points
        // Q_i = p * (P_{i+1} - P_i) / (U_{i+p+1} - U_{i+1})
        // @cite biagiotti2009 -- Sec. 4.5
        auto const vel = de_boor(tc, d1_points_, d1_knots_, p - 1);

        // Acceleration: second derivative, also precomputed at construction
        auto acc = Scalar{0};
        if constexpr (p >= 2) {
            acc = de_boor(tc, d2_points_, d2_knots_, p - 2);
        }

        return {
            .position = Vector<Scalar, 1>{pos},
            .velocity = Vector<Scalar, 1>{vel},
            .acceleration = Vector<Scalar, 1>{acc},
        };
    }

    /// @brief Total trajectory duration (active parameter range).
    ///
    /// Returns U[n+1] - U[p] where n = control_points.size()-1 and p = Degree.
    auto duration() const -> Scalar
    {
        constexpr int p = Degree;
        auto const n = static_cast<int>(control_points_.size()) - 1;
        return knots_[static_cast<std::size_t>(n + 1)] - knots_[static_cast<std::size_t>(p)];
    }

  private:
    /// @brief Tag selecting the non-validating constructor reserved for `create`.
    struct unchecked_t
    {
        explicit unchecked_t() = default;
    };

    /// @brief Construct from a configuration already validated by `create`.
    ///
    /// If knot_vector is empty, generates a uniform clamped knot vector:
    /// m = n + p + 1 total knots, first p+1 = 0, last p+1 = 1, interior uniform.
    ///
    /// @cite biagiotti2009 -- Sec. 4.5
    bspline_trajectory(unchecked_t, config const& cfg)
        : control_points_{cfg.control_points}
    {
        auto const n = static_cast<int>(control_points_.size()) - 1;
        constexpr int p = Degree;

        if (cfg.knot_vector.empty()) {
            // Auto-generate uniform clamped knot vector
            // m + 1 = n + p + 2 total knot values
            // @cite biagiotti2009 -- Sec. 4.5
            auto const m = n + p + 1;
            knots_.resize(static_cast<std::size_t>(m + 1));

            // First p+1 knots = 0
            for (int i = 0; i <= p; ++i) {
                knots_[static_cast<std::size_t>(i)] = Scalar{0};
            }
            // Last p+1 knots = 1
            for (int i = m - p; i <= m; ++i) {
                knots_[static_cast<std::size_t>(i)] = Scalar{1};
            }
            // Interior knots uniformly spaced
            auto const num_interior = n - p;
            for (int i = 1; i <= num_interior; ++i) {
                knots_[static_cast<std::size_t>(p + i)] =
                    static_cast<Scalar>(i) / static_cast<Scalar>(num_interior + 1);
            }
        } else {
            knots_ = cfg.knot_vector;
        }

        // Precompute the derivative control points and knot vectors once so
        // evaluate() stays allocation free on the per-tick hot path. They
        // depend only on control_points_ and knots_, which are fixed here.
        auto d1 = derivative_data(control_points_, knots_, p);
        d1_points_ = std::move(d1.first);
        d1_knots_ = std::move(d1.second);

        if constexpr (p >= 2) {
            auto d2 = derivative_data(d1_points_, d1_knots_, p - 1);
            d2_points_ = std::move(d2.first);
            d2_knots_ = std::move(d2.second);
        }
    }

    std::vector<Scalar> control_points_;
    std::vector<Scalar> knots_;
    std::vector<Scalar> d1_points_;
    std::vector<Scalar> d1_knots_;
    std::vector<Scalar> d2_points_;
    std::vector<Scalar> d2_knots_;

    /// @brief De Boor's algorithm for B-spline evaluation.
    ///
    /// Finds span k such that U[k] <= t < U[k+1], then performs triangular
    /// computation to produce the point on the curve.
    ///
    /// @cite biagiotti2009 -- Sec. 4.5
    static auto de_boor(
        Scalar t,
        std::vector<Scalar> const& points,
        std::vector<Scalar> const& U,
        int degree) -> Scalar
    {
        auto const n = static_cast<int>(points.size()) - 1;

        // Find span k such that U[k] <= t < U[k+1]
        // For the last knot, we use k = n (clamped endpoint)
        int k = find_span(t, U, n, degree);

        // Initialize d[j] = P[j] for j in [k-degree, k]. The buffer is sized
        // by the class degree because every call uses degree <= Degree, so a
        // fixed-size array avoids heap allocation on the evaluation hot path.
        std::array<Scalar, static_cast<std::size_t>(Degree) + 1> d{};
        for (int j = 0; j <= degree; ++j) {
            d[static_cast<std::size_t>(j)] = points[static_cast<std::size_t>(k - degree + j)];
        }

        // Triangular computation
        for (int r = 1; r <= degree; ++r) {
            for (int j = degree; j >= r; --j) {
                auto const idx = k - degree + j;
                auto const denom =
                    U[static_cast<std::size_t>(idx + degree - r + 1)] -
                    U[static_cast<std::size_t>(idx)];
                if (std::abs(denom) < Scalar{1e-10} * (Scalar{1} + std::abs(U[static_cast<std::size_t>(idx)]))) {
                    // Near-zero knot span relative to knot magnitude, keep value
                    continue;
                }
                auto const alpha = (t - U[static_cast<std::size_t>(idx)]) / denom;
                d[static_cast<std::size_t>(j)] =
                    (Scalar{1} - alpha) * d[static_cast<std::size_t>(j - 1)] +
                    alpha * d[static_cast<std::size_t>(j)];
            }
        }

        return d[static_cast<std::size_t>(degree)];
    }

    /// @brief Find the knot span index k such that U[k] <= t < U[k+1].
    static auto find_span(
        Scalar t,
        std::vector<Scalar> const& U,
        int n,
        int degree) -> int
    {
        // Handle endpoint: t == U[n+1]
        if (t >= U[static_cast<std::size_t>(n + 1)]) {
            return n;
        }

        // Binary search via upper_bound
        auto it = std::upper_bound(
            U.begin() + degree,
            U.begin() + n + 2,
            t);

        auto k = static_cast<int>(it - U.begin()) - 1;

        // Clamp to valid range [degree, n]
        k = std::clamp(k, degree, n);
        return k;
    }

    /// @brief Compute derivative control points and knot vector.
    ///
    /// Q_i = deg * (P_{i+1} - P_i) / (U_{i+deg+1} - U_{i+1})
    /// The derivative knot vector drops the first and last knot.
    ///
    /// @cite biagiotti2009 -- Sec. 4.5
    static auto derivative_data(
        std::vector<Scalar> const& points,
        std::vector<Scalar> const& U,
        int deg) -> std::pair<std::vector<Scalar>, std::vector<Scalar>>
    {
        auto const n = static_cast<int>(points.size()) - 1;

        std::vector<Scalar> d_points(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            auto const denom =
                U[static_cast<std::size_t>(i + deg + 1)] -
                U[static_cast<std::size_t>(i + 1)];
            if (std::abs(denom) < Scalar{1e-10} * (Scalar{1} + std::abs(U[static_cast<std::size_t>(i + 1)]))) {
                d_points[static_cast<std::size_t>(i)] = Scalar{0};
            } else {
                d_points[static_cast<std::size_t>(i)] =
                    static_cast<Scalar>(deg) *
                    (points[static_cast<std::size_t>(i + 1)] -
                     points[static_cast<std::size_t>(i)]) /
                    denom;
            }
        }

        // Derivative knot vector: drop first and last knot
        std::vector<Scalar> d_knots(U.begin() + 1, U.end() - 1);

        return {d_points, d_knots};
    }
};

namespace detail
{

/// @brief Evaluate B-spline basis function B_{i,p}(t) using Cox-de Boor recursion.
///
/// @cite biagiotti2009 -- Sec. 4.5
template <typename Scalar>
auto bspline_basis(
    int i,
    int p,
    Scalar t,
    std::vector<Scalar> const& U) -> Scalar
{
    if (p == 0) {
        auto const ui = U[static_cast<std::size_t>(i)];
        auto const ui1 = U[static_cast<std::size_t>(i + 1)];
        // Standard: B_{i,0}(t) = 1 if U[i] <= t < U[i+1]
        // At the right endpoint of the entire domain, include the right
        // boundary for any span with non-zero width (clamped endpoint rule).
        auto const t_max = U.back();
        if (t >= ui && t < ui1) {
            return Scalar{1};
        }
        if (t == t_max && ui < ui1 && t >= ui && t <= ui1) {
            return Scalar{1};
        }
        return Scalar{0};
    }

    Scalar left{0};
    Scalar right{0};

    auto const denom_left =
        U[static_cast<std::size_t>(i + p)] - U[static_cast<std::size_t>(i)];
    if (std::abs(denom_left) > std::numeric_limits<Scalar>::epsilon()) {
        left = (t - U[static_cast<std::size_t>(i)]) / denom_left *
               bspline_basis(i, p - 1, t, U);
    }

    auto const denom_right =
        U[static_cast<std::size_t>(i + p + 1)] - U[static_cast<std::size_t>(i + 1)];
    if (std::abs(denom_right) > std::numeric_limits<Scalar>::epsilon()) {
        right = (U[static_cast<std::size_t>(i + p + 1)] - t) / denom_right *
                bspline_basis(i + 1, p - 1, t, U);
    }

    return left + right;
}

}

/// @brief Evaluate B-spline basis function B_{i,p}(t) (free function for interpolation).
///
/// @cite biagiotti2009 -- Sec. 4.5
template <typename Scalar>
auto basis_function(
    int i,
    int p,
    Scalar t,
    std::vector<Scalar> const& U) -> Scalar
{
    return detail::bspline_basis(i, p, t, U);
}

/// @brief Factory: construct B-spline that interpolates waypoints at given parameter values.
///
/// Generates a clamped knot vector from the parameter values, builds the
/// interpolation matrix of basis function values N[i][j], and solves N * P = Q
/// for the control points P.
///
/// Rejections, checked in order:
///  * times/positions length mismatch      -> spline_error::size_mismatch
///  * fewer than Degree + 1 waypoints      -> spline_error::too_few_points
/// Any downstream `bspline_trajectory::create` failure is propagated.
///
/// @tparam Scalar  Floating-point type
/// @tparam Degree  B-spline degree (compile-time)
/// @param times     n+1 parameter values (must be strictly increasing)
/// @param positions n+1 waypoint positions
/// @return expected carrying a bspline_trajectory that passes through all waypoints
///
/// @cite biagiotti2009 -- Sec. 4.5
template <typename Scalar, int Degree>
auto make_bspline_interpolation(
    std::vector<Scalar> const& times,
    std::vector<Scalar> const& positions)
    -> ctrlpp::expected<bspline_trajectory<Scalar, Degree>, spline_error>
{
    if (times.size() != positions.size()) {
        return ctrlpp::unexpected(spline_error::size_mismatch);
    }
    if (times.size() < static_cast<std::size_t>(Degree + 1)) {
        return ctrlpp::unexpected(spline_error::too_few_points);
    }

    auto const n = static_cast<int>(times.size()) - 1;
    constexpr int p = Degree;

    // Generate clamped knot vector from parameter values
    // First p+1 knots = times[0], last p+1 knots = times[n]
    // Interior knots: averaging per de Boor's method
    // u_{j+p} = (1/p) * sum_{i=j}^{j+p-1} t_i  for j = 1, ..., n-p
    auto const m = n + p + 1;
    std::vector<Scalar> knots(static_cast<std::size_t>(m + 1));

    for (int i = 0; i <= p; ++i) {
        knots[static_cast<std::size_t>(i)] = times[0];
    }
    for (int i = m - p; i <= m; ++i) {
        knots[static_cast<std::size_t>(i)] = times[static_cast<std::size_t>(n)];
    }

    // Interior knots via averaging
    for (int j = 1; j <= n - p; ++j) {
        Scalar sum{0};
        for (int i = j; i < j + p; ++i) {
            sum += times[static_cast<std::size_t>(i)];
        }
        knots[static_cast<std::size_t>(j + p)] = sum / static_cast<Scalar>(p);
    }

    // Build interpolation matrix N[i][j] = B_{j,p}(times[i])
    // using de Boor-Cox recursive basis function evaluation
    auto const num_basis = n + 1;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> N(num_basis, num_basis);
    N.setZero();

    for (int i = 0; i <= n; ++i) {
        auto const t = times[static_cast<std::size_t>(i)];
        for (int j = 0; j <= n; ++j) {
            N(i, j) = basis_function(j, p, t, knots);
        }
    }

    // Right-hand side
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> rhs(num_basis);
    for (int i = 0; i <= n; ++i) {
        rhs(i) = positions[static_cast<std::size_t>(i)];
    }

    // Solve for control points
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> ctrl = N.colPivHouseholderQr().solve(rhs);

    std::vector<Scalar> control_points(static_cast<std::size_t>(num_basis));
    for (int i = 0; i < num_basis; ++i) {
        control_points[static_cast<std::size_t>(i)] = ctrl(i);
    }

    return bspline_trajectory<Scalar, Degree>::create({
        .control_points = std::move(control_points),
        .knot_vector = std::move(knots),
    });
}

static_assert(trajectory_segment<bspline_trajectory<double, 3>, double, 1>);
static_assert(trajectory_segment<bspline_trajectory<double, 5>, double, 1>);

}

#endif
