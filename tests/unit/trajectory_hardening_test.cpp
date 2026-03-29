#include "hardening_helpers.h"

#include "ctrlpp/trajectory/cubic_spline.h"
#include "ctrlpp/trajectory/smoothing_spline.h"
#include "ctrlpp/trajectory/bspline_trajectory.h"
#include "ctrlpp/trajectory/trapezoidal_trajectory.h"
#include "ctrlpp/trajectory/double_s_trajectory.h"
#include "ctrlpp/trajectory/online_planner_2nd.h"
#include "ctrlpp/trajectory/online_planner_3rd.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

using Catch::Matchers::WithinAbs;

// ── Cubic spline hardening ─────────────────────────────────────────────────────

TEST_CASE("Cubic spline with exactly 2 points", "[cubic_spline][hardening][negative]")
{
    ctrlpp::cubic_spline<double>::config cfg{
        .times = {0.0, 1.0},
        .positions = {0.0, 1.0},
    };

    ctrlpp::cubic_spline<double> spline(cfg);
    auto pt = spline.evaluate(0.5);
    REQUIRE(std::isfinite(pt.position(0)));
    REQUIRE(std::isfinite(pt.velocity(0)));
}

TEST_CASE("Cubic spline interpolation matches at knots", "[cubic_spline][hardening][precision]")
{
    std::vector<double> times{0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> positions{0.0, 1.0, 0.5, 2.0, 1.5};

    ctrlpp::cubic_spline<double>::config cfg{
        .times = times,
        .positions = positions,
    };

    ctrlpp::cubic_spline<double> spline(cfg);

    for (std::size_t i = 0; i < times.size(); ++i) {
        auto pt = spline.evaluate(times[i]);
        REQUIRE_THAT(pt.position(0), WithinAbs(positions[i], 1e-10));
    }
}

TEST_CASE("Cubic spline with huge span", "[cubic_spline][hardening][negative]")
{
    ctrlpp::cubic_spline<double>::config cfg{
        .times = {0.0, 1e15, 2e15},
        .positions = {0.0, 1.0, 0.0},
    };

    ctrlpp::cubic_spline<double> spline(cfg);
    auto pt = spline.evaluate(0.5e15);
    REQUIRE(std::isfinite(pt.position(0)));
}

// ── Smoothing spline hardening ─────────────────────────────────────────────────

TEST_CASE("Smoothing spline with mu=1 approaches interpolation", "[smoothing_spline][hardening][negative]")
{
    std::vector<double> times{0.0, 1.0, 2.0, 3.0};
    std::vector<double> positions{0.0, 1.0, 0.5, 2.0};

    ctrlpp::smoothing_spline<double>::config cfg{
        .times = times,
        .positions = positions,
        .mu = 1.0,
    };

    ctrlpp::smoothing_spline<double> spline(cfg);

    // At mu=1, should pass through all waypoints
    for (std::size_t i = 0; i < times.size(); ++i) {
        auto pt = spline.evaluate(times[i]);
        REQUIRE_THAT(pt.position(0), WithinAbs(positions[i], 0.01));
    }
}

TEST_CASE("Smoothing spline with very large lambda", "[smoothing_spline][hardening][negative]")
{
    std::vector<double> times{0.0, 1.0, 2.0, 3.0};
    std::vector<double> positions{0.0, 1.0, 0.5, 2.0};

    // mu near zero -> lambda very large -> maximum smoothness (near straight line)
    ctrlpp::smoothing_spline<double>::config cfg{
        .times = times,
        .positions = positions,
        .mu = 0.001,
    };

    ctrlpp::smoothing_spline<double> spline(cfg);
    auto pt = spline.evaluate(1.5);
    REQUIRE(std::isfinite(pt.position(0)));
}

// ── B-spline hardening ─────────────────────────────────────────────────────────

TEST_CASE("B-spline with insufficient control points", "[bspline][hardening][negative]")
{
    // Degree 3 needs at least 4 control points
    ctrlpp::bspline_trajectory<double, 3>::config cfg{
        .control_points = {0.0, 1.0, 2.0}, // Only 3
    };

    REQUIRE_THROWS_AS(ctrlpp::bspline_trajectory<double, 3>(cfg), std::invalid_argument);
}

TEST_CASE("B-spline with non-ascending knot vector", "[bspline][hardening][negative]")
{
    ctrlpp::bspline_trajectory<double, 3>::config cfg{
        .control_points = {0.0, 1.0, 2.0, 3.0, 4.0},
        .knot_vector = {0.0, 0.0, 0.0, 0.0, 0.5, 0.3, 1.0, 1.0, 1.0}, // Non-ascending
    };

    REQUIRE_THROWS(ctrlpp::bspline_trajectory<double, 3>(cfg));
}

// ── Trapezoidal trajectory hardening ───────────────────────────────────────────

TEST_CASE("Trapezoidal with zero distance", "[trapezoidal][hardening][negative]")
{
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 5.0, .q1 = 5.0, .v_max = 1.0, .a_max = 1.0,
    };

    ctrlpp::trapezoidal_trajectory<double> traj(cfg);
    REQUIRE_THAT(traj.duration(), WithinAbs(0.0, 1e-12));
    auto pt = traj.evaluate(0.0);
    REQUIRE_THAT(pt.position(0), WithinAbs(5.0, 1e-12));
}

TEST_CASE("Trapezoidal with negative max velocity", "[trapezoidal][hardening][negative]")
{
    // Negative v_max is unusual but should handle gracefully
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 1.0, .v_max = -1.0, .a_max = 1.0,
    };

    ctrlpp::trapezoidal_trajectory<double> traj(cfg);
    // Should produce finite results (negative v_max may produce degenerate profile)
    REQUIRE(std::isfinite(traj.duration()));
}

TEST_CASE("Trapezoidal triangle profile reaches correct peak velocity", "[trapezoidal][hardening][precision]")
{
    // Short distance forces triangular profile
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 0.5, .v_max = 10.0, .a_max = 2.0,
    };

    ctrlpp::trapezoidal_trajectory<double> traj(cfg);
    REQUIRE(traj.is_triangular());

    // Peak velocity for triangle: sqrt(a * h) = sqrt(2 * 0.5) = 1.0
    REQUIRE_THAT(traj.peak_velocity(), WithinAbs(1.0, 0.01));

    // Final position should be q1
    auto pt = traj.evaluate(traj.duration());
    REQUIRE_THAT(pt.position(0), WithinAbs(0.5, 1e-6));
}

// ── Double-S trajectory hardening ──────────────────────────────────────────────

TEST_CASE("Double-S with zero distance", "[double_s][hardening][negative]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 3.0, .q1 = 3.0, .v_max = 1.0, .a_max = 1.0, .j_max = 1.0,
    };

    ctrlpp::double_s_trajectory<double> traj(cfg);
    REQUIRE_THAT(traj.duration(), WithinAbs(0.0, 1e-12));
    auto pt = traj.evaluate(0.0);
    REQUIRE_THAT(pt.position(0), WithinAbs(3.0, 1e-12));
}

TEST_CASE("Double-S with negative jerk limit", "[double_s][hardening][negative]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 1.0, .v_max = 1.0, .a_max = 1.0, .j_max = -1.0,
    };

    ctrlpp::double_s_trajectory<double> traj(cfg);
    REQUIRE(std::isfinite(traj.duration()));
}

// ── Online planner 2nd hardening ───────────────────────────────────────────────

TEST_CASE("Online planner 2nd with zero max velocity", "[online_planner_2nd][hardening][negative]")
{
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 0.0, .a_max = 1.0};
    ctrlpp::online_planner_2nd<double> planner(cfg);

    planner.update(1.0);
    auto pt = planner.sample(0.1);
    REQUIRE(std::isfinite(pt.position(0)));
}

TEST_CASE("Online planner 2nd with zero max acceleration", "[online_planner_2nd][hardening][negative]")
{
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 1.0, .a_max = 0.0};
    ctrlpp::online_planner_2nd<double> planner(cfg);

    planner.update(1.0);
    auto pt = planner.sample(0.1);
    REQUIRE(std::isfinite(pt.position(0)));
}

TEST_CASE("Online planner 2nd with instant target flip", "[online_planner_2nd][hardening][negative]")
{
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 1.0, .a_max = 2.0};
    ctrlpp::online_planner_2nd<double> planner(cfg);

    planner.update(5.0);
    planner.sample(0.1);
    planner.sample(0.2);

    // Flip target mid-motion
    planner.update(-5.0);
    auto pt = planner.sample(0.3);
    REQUIRE(std::isfinite(pt.position(0)));
    REQUIRE(std::isfinite(pt.velocity(0)));
}

TEST_CASE("Online planner 2nd reaches target", "[online_planner_2nd][hardening][convergence]")
{
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 2.0, .a_max = 1.0};
    ctrlpp::online_planner_2nd<double> planner(cfg);

    planner.update(3.0);

    double t = 0.0;
    for (int i = 0; i < 1000; ++i) {
        t += 0.01;
        planner.sample(t);
    }

    auto pt = planner.sample(t + 0.01);
    REQUIRE_THAT(pt.position(0), WithinAbs(3.0, 1e-6));
}

// ── Online planner 3rd hardening ───────────────────────────────────────────────

TEST_CASE("Online planner 3rd with zero max jerk", "[online_planner_3rd][hardening][negative]")
{
    ctrlpp::online_planner_3rd<double>::config cfg{.v_max = 1.0, .a_max = 1.0, .j_max = 0.0};
    ctrlpp::online_planner_3rd<double> planner(cfg);

    planner.update(1.0);
    auto pt = planner.sample(0.1);
    REQUIRE(std::isfinite(pt.position(0)));
}

TEST_CASE("Online planner 3rd with instant target reversal", "[online_planner_3rd][hardening][negative]")
{
    ctrlpp::online_planner_3rd<double>::config cfg{.v_max = 1.0, .a_max = 2.0, .j_max = 5.0};
    ctrlpp::online_planner_3rd<double> planner(cfg);

    planner.update(5.0);
    planner.sample(0.1);

    // Reverse direction
    planner.update(-5.0);
    auto pt = planner.sample(0.2);
    REQUIRE(std::isfinite(pt.position(0)));
}

TEST_CASE("Online planner 3rd reaches target", "[online_planner_3rd][hardening][convergence]")
{
    ctrlpp::online_planner_3rd<double>::config cfg{.v_max = 2.0, .a_max = 1.0, .j_max = 5.0};
    ctrlpp::online_planner_3rd<double> planner(cfg);

    planner.update(3.0);

    double t = 0.0;
    for (int i = 0; i < 2000; ++i) {
        t += 0.01;
        planner.sample(t);
    }

    auto pt = planner.sample(t + 0.01);
    REQUIRE_THAT(pt.position(0), WithinAbs(3.0, 1e-4));
}
