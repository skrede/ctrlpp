#include <ctrlpp/trajectory/smoothing_spline.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <vector>

using Catch::Matchers::WithinAbs;

TEST_CASE("smoothing_spline mu=1 interpolates exactly", "[smoothing_spline]")
{
    // 5 clean waypoints -- mu=1 should pass through all of them
    std::vector<double> times     = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> positions = {0.0, 1.0, 0.5, 2.0, 1.5};

    ctrlpp::smoothing_spline<double> spline({
        .times = times,
        .positions = positions,
        .mu = 1.0,
    });

    for (std::size_t i = 0; i < times.size(); ++i) {
        auto const pt = spline.evaluate(times[i]);
        REQUIRE_THAT(pt.position[0], WithinAbs(positions[i], 1e-8));
    }
}

TEST_CASE("smoothing_spline mu=0.5 approximates noisy data", "[smoothing_spline]")
{
    // 10 noisy waypoints -- mu=0.5 should NOT pass exactly through them
    // but should stay within noise amplitude
    double const noise_amplitude = 0.35;
    std::vector<double> times     = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5};
    // Underlying signal is sin(t), with added noise
    std::vector<double> positions = {
        0.0 + 0.1, 0.479 - 0.15, 0.841 + 0.2, 0.997 - 0.1, 0.909 + 0.05,
        0.598 - 0.2, 0.141 + 0.15, -0.351 - 0.1, -0.757 + 0.25, -0.978 - 0.05,
    };

    ctrlpp::smoothing_spline<double> spline({
        .times = times,
        .positions = positions,
        .mu = 0.5,
    });

    // Should NOT pass exactly through waypoints
    bool any_deviation = false;
    double max_deviation = 0.0;
    for (std::size_t i = 0; i < times.size(); ++i) {
        auto const pt = spline.evaluate(times[i]);
        double const dev = std::abs(pt.position[0] - positions[i]);
        if (dev > 1e-8) {
            any_deviation = true;
        }
        max_deviation = std::max(max_deviation, dev);
    }
    REQUIRE(any_deviation);
    // Deviation should be bounded by noise amplitude
    REQUIRE(max_deviation < noise_amplitude);
}

TEST_CASE("smoothing_spline mu=0.01 is very smooth", "[smoothing_spline]")
{
    // With very low mu, the spline should deviate significantly from noisy waypoints
    std::vector<double> times     = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> positions = {0.0, 1.5, -0.5, 2.0, -1.0, 0.5};

    ctrlpp::smoothing_spline<double> spline({
        .times = times,
        .positions = positions,
        .mu = 0.01,
    });

    // Should deviate significantly
    double max_dev = 0.0;
    for (std::size_t i = 0; i < times.size(); ++i) {
        auto const pt = spline.evaluate(times[i]);
        max_dev = std::max(max_dev, std::abs(pt.position[0] - positions[i]));
    }
    REQUIRE(max_dev > 0.1);
}

TEST_CASE("smoothing_spline smoothness ordering with mu", "[smoothing_spline]")
{
    // Sum of squared second derivatives should decrease as mu decreases
    std::vector<double> times     = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> positions = {0.0, 2.0, -1.0, 3.0, 0.5};

    auto roughness = [&](double mu) {
        ctrlpp::smoothing_spline<double> spline({
            .times = times,
            .positions = positions,
            .mu = mu,
        });

        double sum_sq_acc = 0.0;
        int const n_samples = 200;
        double const dt = spline.duration() / n_samples;
        for (int k = 0; k <= n_samples; ++k) {
            auto const t = times.front() + k * dt;
            auto const pt = spline.evaluate(t);
            sum_sq_acc += pt.acceleration[0] * pt.acceleration[0];
        }
        return sum_sq_acc * dt;
    };

    double const r_high = roughness(1.0);
    double const r_mid  = roughness(0.5);
    double const r_low  = roughness(0.1);

    REQUIRE(r_high > r_mid);
    REQUIRE(r_mid > r_low);
}

TEST_CASE("smoothing_spline duration", "[smoothing_spline]")
{
    std::vector<double> times     = {1.0, 2.5, 4.0, 6.0};
    std::vector<double> positions = {0.0, 1.0, 0.5, 2.0};

    ctrlpp::smoothing_spline<double> spline({
        .times = times,
        .positions = positions,
        .mu = 0.8,
    });

    REQUIRE_THAT(spline.duration(), WithinAbs(5.0, 1e-12));
}

TEST_CASE("smoothing_spline satisfies trajectory_segment concept", "[smoothing_spline]")
{
    // This is verified by the static_assert in the header,
    // but test it explicitly too
    STATIC_REQUIRE(ctrlpp::trajectory_segment<
        ctrlpp::smoothing_spline<double>, double, 1>);
}

// -- Validation ---------------------------------------------------------------
//
// Before the try_create conversion these invalid configurations were guarded
// only by assert(), so release builds accepted them silently; mu = 0 in
// particular was clamped to an arbitrary small constant instead of rejected.
// Each case pins the specific rejection enumerator.

TEST_CASE("smoothing_spline try_create rejects mu = 0", "[smoothing_spline][validation]")
{
    auto const result = ctrlpp::smoothing_spline<double>::try_create({
        .times = {0.0, 1.0, 2.0},
        .positions = {0.0, 1.0, 0.0},
        .mu = 0.0,
    });

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::mu_out_of_range);
}

TEST_CASE("smoothing_spline try_create rejects mu > 1", "[smoothing_spline][validation]")
{
    auto const result = ctrlpp::smoothing_spline<double>::try_create({
        .times = {0.0, 1.0, 2.0},
        .positions = {0.0, 1.0, 0.0},
        .mu = 2.0,
    });

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::mu_out_of_range);
}

TEST_CASE("smoothing_spline try_create rejects NaN mu", "[smoothing_spline][validation]")
{
    auto const result = ctrlpp::smoothing_spline<double>::try_create({
        .times = {0.0, 1.0, 2.0},
        .positions = {0.0, 1.0, 0.0},
        .mu = std::numeric_limits<double>::quiet_NaN(),
    });

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::mu_out_of_range);
}

TEST_CASE("smoothing_spline try_create accepts mu = 1", "[smoothing_spline][validation]")
{
    std::vector<double> const times     = {0.0, 1.0, 2.0};
    std::vector<double> const positions = {0.0, 1.0, 0.0};

    auto const result = ctrlpp::smoothing_spline<double>::try_create({
        .times = times,
        .positions = positions,
        .mu = 1.0,
    });

    REQUIRE(result.has_value());
    for (std::size_t i = 0; i < times.size(); ++i) {
        REQUIRE_THAT(result->evaluate(times[i]).position[0], WithinAbs(positions[i], 1e-8));
    }
}

TEST_CASE("smoothing_spline try_create rejects fewer than 2 waypoints", "[smoothing_spline][validation]")
{
    auto const result = ctrlpp::smoothing_spline<double>::try_create({
        .times = {0.0},
        .positions = {1.0},
        .mu = 0.5,
    });

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::too_few_points);
}

TEST_CASE("smoothing_spline try_create rejects length mismatch", "[smoothing_spline][validation]")
{
    auto const result = ctrlpp::smoothing_spline<double>::try_create({
        .times = {0.0, 1.0, 2.0},
        .positions = {0.0, 1.0},
        .mu = 0.5,
    });

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::size_mismatch);
}

TEST_CASE("smoothing_spline try_create rejects non-increasing times", "[smoothing_spline][validation]")
{
    auto const result = ctrlpp::smoothing_spline<double>::try_create({
        .times = {0.0, 1.0, 1.0},
        .positions = {0.0, 1.0, 0.0},
        .mu = 0.5,
    });

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::non_increasing_times);
}
