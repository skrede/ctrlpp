#include "ctrlpp/trajectory/bspline_trajectory.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <algorithm>
#include <type_traits>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace
{

// Construction is fallible, so every valid-input site goes through the factory
// and asserts success here; a config that becomes unrealizable fails the test
// instead of quietly skipping it. Rejection cases never use this helper.
template <int Degree>
void require_finite_evaluation(ctrlpp::bspline_trajectory<double, Degree> const& spline, double t_begin)
{
    for (int i = 0; i <= 100; ++i) {
        auto const pt = spline.evaluate(t_begin + spline.duration() * static_cast<double>(i) / 100.0);
        REQUIRE(std::isfinite(pt.position(0)));
        REQUIRE(std::isfinite(pt.velocity(0)));
        REQUIRE(std::isfinite(pt.acceleration(0)));
    }
}

template <int Degree>
auto make_bspline(typename ctrlpp::bspline_trajectory<double, Degree>::config const& cfg)
    -> ctrlpp::bspline_trajectory<double, Degree>
{
    auto created = ctrlpp::bspline_trajectory<double, Degree>::create(cfg);
    REQUIRE(created.has_value());
    auto const t_begin = cfg.knot_vector.empty() ? 0.0 : cfg.knot_vector[static_cast<std::size_t>(Degree)];
    require_finite_evaluation(*created, t_begin);
    return *std::move(created);
}

// Marsden's identity (de Boor 2001, Ch. IX): 3 t^2 + 3 t = sum_i (a b + a c + b c + a + b + c) N_{i,3}(t)
// with (a, b, c) = (U[i + 1], U[i + 2], U[i + 3]).
auto marsden_control_points(std::vector<double> const& knots) -> std::vector<double>
{
    std::vector<double> points(knots.size() - 4);
    for (std::size_t i = 0; i < points.size(); ++i) {
        auto const a = knots[i + 1];
        auto const b = knots[i + 2];
        auto const c = knots[i + 3];
        points[i] = a * b + a * c + b * c + a + b + c;
    }
    return points;
}

}

TEST_CASE("Cubic B-spline endpoint interpolation", "[bspline]")
{
    // Degree 3 with 5 control points and auto-generated uniform clamped knots.
    // Clamped B-spline passes through first and last control points.
    std::vector<double> ctrl = {0.0, 1.0, 3.0, 2.0, 4.0};

    auto bs = make_bspline<3>({.control_points = ctrl});

    auto const p0 = bs.evaluate(0.0);
    auto const pT = bs.evaluate(bs.duration());

    REQUIRE_THAT(p0.position(0), WithinAbs(ctrl.front(), 1e-12));
    REQUIRE_THAT(pT.position(0), WithinAbs(ctrl.back(), 1e-12));
}

TEST_CASE("Cubic B-spline convex hull", "[bspline]")
{
    std::vector<double> ctrl = {1.0, 2.0, 4.0, 3.0, 5.0};

    auto bs = make_bspline<3>({.control_points = ctrl});

    auto const min_cp = *std::min_element(ctrl.begin(), ctrl.end());
    auto const max_cp = *std::max_element(ctrl.begin(), ctrl.end());

    auto const T = bs.duration();
    int const N = 50;
    for (int i = 0; i <= N; ++i) {
        auto const t = T * static_cast<double>(i) / N;
        auto const pt = bs.evaluate(t);
        REQUIRE(pt.position(0) >= min_cp - 1e-12);
        REQUIRE(pt.position(0) <= max_cp + 1e-12);
    }
}

TEST_CASE("make_bspline_interpolation passes through waypoints", "[bspline]")
{
    // 5 waypoints at uniform times
    std::vector<double> times = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> positions = {0.0, 2.0, 1.0, 3.0, 5.0};

    auto const bs = ctrlpp::make_bspline_interpolation<double, 3>(times, positions);
    REQUIRE(bs.has_value());
    require_finite_evaluation(*bs, times.front());

    for (std::size_t i = 0; i < times.size(); ++i) {
        auto const pt = bs->evaluate(times[i]);
        REQUIRE_THAT(pt.position(0), WithinAbs(positions[i], 1e-10));
    }
}

TEST_CASE("Quintic B-spline evaluation", "[bspline]")
{
    // Degree 5 with 7 control points
    std::vector<double> ctrl = {0.0, 1.0, 2.5, 4.0, 3.0, 2.0, 5.0};

    auto bs = make_bspline<5>({.control_points = ctrl});

    auto const T = bs.duration();
    REQUIRE(T > 0.0);

    int const N = 20;
    for (int i = 0; i <= N; ++i) {
        auto const t = T * static_cast<double>(i) / N;
        auto const pt = bs.evaluate(t);
        // Should not produce NaN or Inf
        REQUIRE(std::isfinite(pt.position(0)));
        REQUIRE(std::isfinite(pt.velocity(0)));
        REQUIRE(std::isfinite(pt.acceleration(0)));
    }
}

TEST_CASE("Custom knot vector", "[bspline]")
{
    // Cubic with 5 control points, custom non-uniform knots
    // m = n + p + 1 = 4 + 3 + 1 = 8 knots total (9 values)
    std::vector<double> ctrl = {0.0, 1.0, 3.0, 2.0, 4.0};
    std::vector<double> knots = {0, 0, 0, 0, 0.3, 1, 1, 1, 1};

    auto bs = make_bspline<3>({.control_points = ctrl, .knot_vector = knots});

    auto const p0 = bs.evaluate(0.0);
    auto const pT = bs.evaluate(bs.duration());

    REQUIRE_THAT(p0.position(0), WithinAbs(ctrl.front(), 1e-12));
    REQUIRE_THAT(pT.position(0), WithinAbs(ctrl.back(), 1e-12));
}

// -- Validation ---------------------------------------------------------------
//
// Before the create conversion the constructor threw std::invalid_argument
// for these inputs and make_bspline_interpolation guarded them only by
// assert(), so release builds accepted them silently. Each case pins the
// specific rejection enumerator.

TEST_CASE("create rejects too few control points", "[bspline][validation]")
{
    // Degree 3 needs at least 4 control points
    auto const result = ctrlpp::bspline_trajectory<double, 3>::create({
        .control_points = {0.0, 1.0, 2.0},
    });

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::too_few_control_points);
}

TEST_CASE("create rejects wrong knot vector size", "[bspline][validation]")
{
    // Wrong number of knots for degree 3 with 5 control points (need 9, provide 7)
    auto const result = ctrlpp::bspline_trajectory<double, 3>::create({
        .control_points = {0.0, 1.0, 3.0, 2.0, 4.0},
        .knot_vector = {0, 0, 0, 0.5, 1, 1, 1},
    });

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::bad_knot_count);
}

TEST_CASE("create rejects non-monotonic knot vector", "[bspline][validation]")
{
    // Correct knot count (9) but decreasing interior knots
    auto const result = ctrlpp::bspline_trajectory<double, 3>::create({
        .control_points = {0.0, 1.0, 3.0, 2.0, 4.0},
        .knot_vector = {0, 0, 0, 0, 0.5, 0.3, 1, 1, 1},
    });

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::non_monotonic_knots);
}

TEST_CASE("make_bspline_interpolation rejects length mismatch", "[bspline][validation]")
{
    std::vector<double> times = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> positions = {0.0, 2.0, 1.0, 3.0};

    auto const result = ctrlpp::make_bspline_interpolation<double, 3>(times, positions);

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::size_mismatch);
}

TEST_CASE("make_bspline_interpolation rejects too few waypoints", "[bspline][validation]")
{
    // Degree 3 interpolation needs at least 4 waypoints
    std::vector<double> times = {0.0, 1.0, 2.0};
    std::vector<double> positions = {0.0, 2.0, 1.0};

    auto const result = ctrlpp::make_bspline_interpolation<double, 3>(times, positions);

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::too_few_points);
}

TEST_CASE("Duration returns active parameter range", "[bspline]")
{
    // For uniform clamped knots on [0,1], duration should be 1.0
    std::vector<double> ctrl = {0.0, 1.0, 3.0, 2.0, 4.0};

    auto bs = make_bspline<3>({.control_points = ctrl});

    REQUIRE_THAT(bs.duration(), WithinAbs(1.0, 1e-12));
}

TEST_CASE("Satisfies trajectory_segment concept", "[bspline]")
{
    static_assert(ctrlpp::trajectory_segment<ctrlpp::bspline_trajectory<double, 3>, double, 1>);
    static_assert(ctrlpp::trajectory_segment<ctrlpp::bspline_trajectory<double, 5>, double, 1>);
    SUCCEED("Static assert passed");
}

TEST_CASE("Velocity and acceleration via derivative de Boor", "[bspline]")
{
    // For a linear B-spline (control points on a line), velocity should be constant
    // and acceleration should be zero.
    std::vector<double> ctrl = {0.0, 1.0, 2.0, 3.0, 4.0};

    auto bs = make_bspline<3>({.control_points = ctrl});

    auto const T = bs.duration();
    int const N = 20;
    for (int i = 1; i < N; ++i) {
        auto const t = T * static_cast<double>(i) / N;
        auto const pt = bs.evaluate(t);
        // For collinear control points with uniform clamped knots,
        // velocity should be approximately constant at interior points
        REQUIRE(std::isfinite(pt.velocity(0)));
        REQUIRE(std::isfinite(pt.acceleration(0)));
    }

    // At interior points, velocity should be positive for ascending control points
    auto const mid = bs.evaluate(T * 0.5);
    REQUIRE(mid.velocity(0) > 0.0);
}

TEST_CASE("make_bspline_interpolation rejects a repeated time", "[trajectory][bspline][validation]")
{
    std::vector<double> times = {0.0, 1.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> positions = {0.0, 1.0, 5.0, 2.0, 3.0, 4.0};

    auto const result = ctrlpp::make_bspline_interpolation<double, 3>(times, positions);

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::non_increasing_times);
}

TEST_CASE("make_bspline_interpolation rejects a decreasing time", "[trajectory][bspline][validation]")
{
    std::vector<double> times = {0.0, 2.0, 1.0, 3.0, 4.0};
    std::vector<double> positions = {0.0, 1.0, 2.0, 3.0, 4.0};

    auto const result = ctrlpp::make_bspline_interpolation<double, 3>(times, positions);

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::non_increasing_times);
}

TEST_CASE("make_bspline_interpolation rejects non-finite input", "[trajectory][bspline][validation]")
{
    auto const nan = std::numeric_limits<double>::quiet_NaN();
    auto const inf = std::numeric_limits<double>::infinity();
    std::vector<double> const times = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> const positions = {0.0, 2.0, 1.0, 3.0, 5.0};

    auto with = [](std::vector<double> values, double bad) {
        values[2] = bad;
        return values;
    };

    auto const nan_time = ctrlpp::make_bspline_interpolation<double, 3>(with(times, nan), positions);
    auto const inf_time = ctrlpp::make_bspline_interpolation<double, 3>(with(times, inf), positions);
    auto const nan_position = ctrlpp::make_bspline_interpolation<double, 3>(times, with(positions, nan));

    REQUIRE_FALSE(nan_time.has_value());
    REQUIRE(nan_time.error() == ctrlpp::spline_error::non_finite_input);
    REQUIRE_FALSE(inf_time.has_value());
    REQUIRE(inf_time.error() == ctrlpp::spline_error::non_finite_input);
    REQUIRE_FALSE(nan_position.has_value());
    REQUIRE(nan_position.error() == ctrlpp::spline_error::non_finite_input);
}

TEST_CASE("make_bspline_interpolation rejects a time span that overflows", "[trajectory][bspline][validation]")
{
    std::vector<double> times = {-1.5e308, 1.5e308, 1.6e308, 1.7e308};
    std::vector<double> positions = {0.0, 1.0, 2.0, 3.0};

    auto const result = ctrlpp::make_bspline_interpolation<double, 3>(times, positions);

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::unrepresentable_spline);
}

TEST_CASE("make_bspline_interpolation rejects control points that leave the range", "[trajectory][bspline][validation]")
{
    std::vector<double> times = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> positions = {0.0, 1.0, 1e308, 2.0, 3.0, 4.0};

    auto const result = ctrlpp::make_bspline_interpolation<double, 3>(times, positions);

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::unrepresentable_spline);
}

TEST_CASE("make_bspline_interpolation rejects a time range whose basis overflows", "[trajectory][bspline][validation]")
{
    std::vector<double> times = {-1e308, -0.5e308, 0.0, 0.5e308, 1e308};
    std::vector<double> positions = {0.0, 1.0, 2.0, 3.0, 4.0};

    auto const result = ctrlpp::make_bspline_interpolation<double, 3>(times, positions);

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::unrepresentable_spline);
}

TEST_CASE("create rejects a NaN control point", "[trajectory][bspline][validation]")
{
    auto const result = ctrlpp::bspline_trajectory<double, 3>::create({
        .control_points = {0.0, 1.0, std::numeric_limits<double>::quiet_NaN(), 2.0, 4.0},
    });

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::non_finite_input);
}

TEST_CASE("create rejects infinite knots", "[trajectory][bspline][validation]")
{
    auto const inf = std::numeric_limits<double>::infinity();
    auto const result = ctrlpp::bspline_trajectory<double, 3>::create({
        .control_points = {0.0, 1.0, 2.0, 3.0},
        .knot_vector = {0.0, 0.0, 0.0, 0.0, inf, inf, inf, inf},
    });

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::non_finite_input);
}

TEST_CASE("create rejects control points whose derivative overflows", "[trajectory][bspline][validation]")
{
    auto const result = ctrlpp::bspline_trajectory<double, 3>::create({
        .control_points = {0.0, 1.0, 1.5e308, -1.5e308},
    });

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::unrepresentable_spline);
}

TEMPLATE_TEST_CASE("create rejects a knot vector whose range overflows", "[trajectory][bspline][validation]", float, double)
{
    using Scalar = TestType;
    auto const m = std::numeric_limits<Scalar>::max() / static_cast<Scalar>(1.5);
    auto const result = ctrlpp::bspline_trajectory<Scalar, 3>::create({
        .control_points = {Scalar{0}, Scalar{1}, Scalar{2}, Scalar{3}, Scalar{4}},
        .knot_vector = {-m, -m, -m, -m, Scalar{0}, m, m, m, m},
    });

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::unrepresentable_spline);
}

TEST_CASE("A knot span far below the knot scale is evaluated as given", "[trajectory][bspline]")
{
    // The curve is 3 t^2 + 3 t on [-1, 1]. These dyadic knots make every control point exact, and the
    // central difference of a quadratic is exact for any h, so only rounding separates it from v.
    auto const delta = std::ldexp(1.0, -35);
    std::vector<std::vector<double>> const knot_vectors = {
        {-1, -1, -1, -1, 0, delta, 1, 1, 1, 1},
        {-1, -1, -1, -1, 0, delta, 2 * delta, 3 * delta, 1, 1, 1, 1},
    };
    // A de Boor level adds at most 7 u max|d|: alpha's two roundings weight a difference of at most
    // 2 max|d|, then one each for 1 - alpha, the two products together, and the sum. Position takes
    // 3 levels on |P| <= 6; velocity takes 4 roundings in Q and 2 levels on |Q| <= 9; the difference
    // adds the rounding of q(t + h) - q(t - h), at most 2 max|P|, all divided by the exact 2 h.
    auto const u = std::numeric_limits<double>::epsilon() / 2.0;
    auto const h = 0.0625;
    auto const bound = 44.0 * u * 6.0 / (2.0 * h) + 18.0 * u * 9.0;

    std::vector<double> grid;
    for (std::int32_t j = -15; j <= 15; ++j) {
        grid.push_back(static_cast<double>(j) / 16.0);
    }
    for (std::int32_t k = -64; k <= 160; ++k) {
        grid.push_back(std::ldexp(static_cast<double>(k), -40));
    }

    for (auto const& knots : knot_vectors) {
        auto const spline = make_bspline<3>({.control_points = marsden_control_points(knots), .knot_vector = knots});
        for (auto const t : grid) {
            auto const difference =
                (spline.evaluate(t + h).position(0) - spline.evaluate(t - h).position(0)) / (2.0 * h);
            INFO("knots " << knots.size() << " t " << t);
            REQUIRE_THAT(spline.evaluate(t).velocity(0), WithinAbs(difference, bound));
        }
    }
}

TEMPLATE_TEST_CASE("Knots far below unit scale evaluate the same normalized curve", "[trajectory][bspline]", float, double)
{
    using Scalar = TestType;
    // Scaling every knot by 2^-s commutes with rounding while nothing goes subnormal or overflows, so the
    // scaled curve is a bit-exact image. float interpolates at 2^-40 because the second derivative
    // scales by 2^(2 s) and 2^132 exceeds FLT_MAX.
    std::int32_t const basis_shift = 66;
    std::int32_t const shift = std::is_same_v<Scalar, float> ? 40 : basis_shift;
    auto const scaled = [](std::vector<Scalar> values, std::int32_t s) {
        for (auto& value : values) {
            value = std::ldexp(value, -s);
        }
        return values;
    };

    std::vector<Scalar> const knots = {0, 0, 0, 0, Scalar{0.25F}, Scalar{0.5F}, 1, 1, 1, 1};
    auto const small_knots = scaled(knots, basis_shift);
    for (std::int32_t i = 0; i <= 5; ++i) {
        for (std::int32_t k = 0; k <= 16; ++k) {
            auto const t = static_cast<Scalar>(k) / Scalar{16};
            REQUIRE(ctrlpp::basis_function(i, 3, std::ldexp(t, -basis_shift), small_knots) ==
                    ctrlpp::basis_function(i, 3, t, knots));
        }
    }

    std::vector<Scalar> const times = {0, Scalar{0.25F}, Scalar{0.5F}, Scalar{0.625F}, Scalar{0.75F}, 1};
    std::vector<Scalar> const positions = {0, 1, 3, 2, 4, 5};
    auto const unit = ctrlpp::make_bspline_interpolation<Scalar, 3>(times, positions);
    auto const small = ctrlpp::make_bspline_interpolation<Scalar, 3>(scaled(times, shift), positions);
    REQUIRE(unit.has_value());
    REQUIRE(small.has_value());
    for (std::int32_t k = 0; k <= 64; ++k) {
        auto const t = static_cast<Scalar>(k) / Scalar{64};
        auto const expected = unit->evaluate(t);
        auto const actual = small->evaluate(std::ldexp(t, -shift));
        REQUIRE(actual.position(0) == expected.position(0));
        REQUIRE(actual.velocity(0) == std::ldexp(expected.velocity(0), shift));
        REQUIRE(actual.acceleration(0) == std::ldexp(expected.acceleration(0), 2 * shift));
    }
}

TEST_CASE("An empty knot span at the end of the domain evaluates finite", "[trajectory][bspline]")
{
    // Only t == U[n + 1] with U[n] == U[n + 1] reaches an empty span inside de Boor's triangle.
    std::vector<std::vector<double>> const knot_vectors = {
        {0, 1, 2, 3, 4, 4, 6, 7, 8},
        {2, 2, 2, 2, 2, 2, 2, 2, 2},
    };
    for (auto const& knots : knot_vectors) {
        auto const spline = make_bspline<3>({.control_points = {0.0, 1.0, 3.0, 2.0, 4.0}, .knot_vector = knots});
        auto const end = spline.evaluate(knots[5]);
        REQUIRE(std::isfinite(end.position(0)));
        REQUIRE(std::isfinite(end.velocity(0)));
        REQUIRE(std::isfinite(end.acceleration(0)));
    }
}

TEMPLATE_TEST_CASE("A domain ending on an empty knot span evaluates its left limit", "[trajectory][bspline]", float, double)
{
    using Scalar = TestType;
    std::vector<Scalar> const control_points = {0, 1, 3, 2, 4};

    // B_{4,3} has the empty support [1, 1], so the curve is the Bezier cubic on P0..P3:
    // q(1) = P3, v(1) = 3 (P3 - P2), a(1) = 6 (P3 - 2 P2 + P1), all exact in both types.
    auto const clamped = ctrlpp::bspline_trajectory<Scalar, 3>::create({
        .control_points = control_points,
        .knot_vector = {0, 0, 0, 0, 1, 1, 1, 1, 1},
    });
    REQUIRE(clamped.has_value());
    auto const at_one = clamped->evaluate(Scalar{1});
    REQUIRE(at_one.position(0) == Scalar{2});
    REQUIRE(at_one.velocity(0) == Scalar{-3});
    REQUIRE(at_one.acceleration(0) == Scalar{-18});

    // Second-derivative control points are {2, -8, 3} on knots {2, 3, 4, 4, 6}; the span [3, 4] holds -8.
    auto const double_knot = ctrlpp::bspline_trajectory<Scalar, 3>::create({
        .control_points = control_points,
        .knot_vector = {0, 1, 2, 3, 4, 4, 6, 7, 8},
    });
    REQUIRE(double_knot.has_value());
    auto const at_four = double_knot->evaluate(Scalar{4});
    REQUIRE(at_four.velocity(0) == Scalar{-1});
    REQUIRE(at_four.acceleration(0) == Scalar{-8});
}

TEMPLATE_TEST_CASE("basis_function stays a partition of unity beside a span of the smallest normal width", "[trajectory][bspline]", float, double)
{
    using Scalar = TestType;
    // (t - U[3]) / (U[4] - U[3]) = 8 / min() overflows, where B_{3,0}(t) is zero.
    auto const tiny = std::numeric_limits<Scalar>::min();
    std::vector<Scalar> const knots = {0, 0, 0, 0, tiny, 8, 8, 8, 8};
    for (auto const t : {Scalar{8}, Scalar{4}, tiny}) {
        Scalar sum{0};
        for (std::int32_t i = 0; i <= 4; ++i) {
            auto const value = ctrlpp::basis_function(i, 3, t, knots);
            INFO("i " << i << " t " << t);
            REQUIRE(value >= Scalar{0});
            sum += value;
        }
        REQUIRE_THAT(static_cast<double>(sum), WithinAbs(1.0, 4.0 * static_cast<double>(std::numeric_limits<Scalar>::epsilon())));
    }
}
