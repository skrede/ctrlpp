#include "ctrlpp/trajectory/cubic_spline.h"
#include "ctrlpp/trajectory/detail/tridiagonal.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <vector>
#include <utility>

using Catch::Matchers::WithinAbs;

namespace
{

// Construction is fallible, so every valid-input site goes through the factory
// and asserts success here; a config that becomes unrealizable fails the test
// instead of quietly skipping it. Rejection cases never use this helper.
auto make_spline(ctrlpp::cubic_spline<double>::config const& cfg) -> ctrlpp::cubic_spline<double>
{
    auto created = ctrlpp::cubic_spline<double>::create(cfg);
    REQUIRE(created.has_value());
    return *std::move(created);
}

auto solve_tridiagonal(std::vector<double> const& a, std::vector<double> b, std::vector<double> const& c)
    -> ctrlpp::expected<void, ctrlpp::spline_error>
{
    std::vector<double> d(b.size(), 1.0);
    return ctrlpp::detail::thomas_solve(a, b, c, d);
}

}

// -- Tridiagonal solver tests -------------------------------------------------

TEST_CASE("thomas_solve: 3x3 system", "[trajectory][tridiagonal]")
{
    // System: [[2,1,0],[1,3,1],[0,1,2]] * x = [1,2,3]
    std::vector<double> a = {0.0, 1.0, 1.0};    // sub-diagonal (a[0] unused)
    std::vector<double> b = {2.0, 3.0, 2.0};    // main diagonal
    std::vector<double> c = {1.0, 1.0, 0.0};    // super-diagonal (c[n-1] unused)
    std::vector<double> d = {1.0, 2.0, 3.0};    // RHS

    // Verify against Eigen dense solve
    Eigen::Matrix3d A;
    A << 2, 1, 0,
         1, 3, 1,
         0, 1, 2;
    Eigen::Vector3d rhs(1.0, 2.0, 3.0);
    Eigen::Vector3d x_ref = A.colPivHouseholderQr().solve(rhs);

    auto const solved = ctrlpp::detail::thomas_solve(a, b, c, d);
    REQUIRE(solved.has_value());

    for (int i = 0; i < 3; ++i) {
        REQUIRE_THAT(d[static_cast<std::size_t>(i)], WithinAbs(x_ref(i), 1e-12));
    }
}

TEST_CASE("thomas_solve: 5x5 diagonal-dominant system", "[trajectory][tridiagonal]")
{
    constexpr int n = 5;
    std::vector<double> a(n), b(n), c(n), d(n);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd rhs(n);

    for (int i = 0; i < n; ++i) {
        b[static_cast<std::size_t>(i)] = 4.0;
        A(i, i) = 4.0;
        if (i > 0) {
            a[static_cast<std::size_t>(i)] = 1.0;
            A(i, i - 1) = 1.0;
        }
        if (i < n - 1) {
            c[static_cast<std::size_t>(i)] = 1.0;
            A(i, i + 1) = 1.0;
        }
        d[static_cast<std::size_t>(i)] = static_cast<double>(i + 1);
        rhs(i) = static_cast<double>(i + 1);
    }

    Eigen::VectorXd x_ref = A.colPivHouseholderQr().solve(rhs);

    auto const solved = ctrlpp::detail::thomas_solve(a, b, c, d);
    REQUIRE(solved.has_value());

    for (int i = 0; i < n; ++i) {
        REQUIRE_THAT(d[static_cast<std::size_t>(i)], WithinAbs(x_ref(i), 1e-12));
    }
}

TEST_CASE("cyclic_thomas_solve: 3x3 cyclic system", "[trajectory][tridiagonal]")
{
    constexpr int n = 3;
    std::vector<double> a = {0.0, 1.0, 1.0};
    std::vector<double> b = {4.0, 4.0, 4.0};
    std::vector<double> c = {1.0, 1.0, 0.0};
    std::vector<double> d = {1.0, 2.0, 3.0};
    double alpha = 0.5;  // A[0][n-1]
    double beta = 0.7;   // A[n-1][0]

    Eigen::Matrix3d A;
    A << 4.0, 1.0, 0.5,
         1.0, 4.0, 1.0,
         0.7, 1.0, 4.0;
    Eigen::Vector3d rhs(1.0, 2.0, 3.0);
    Eigen::Vector3d x_ref = A.colPivHouseholderQr().solve(rhs);

    auto const solved = ctrlpp::detail::cyclic_thomas_solve(a, b, c, d, alpha, beta);
    REQUIRE(solved.has_value());

    for (int i = 0; i < n; ++i) {
        REQUIRE_THAT(d[static_cast<std::size_t>(i)], WithinAbs(x_ref(i), 1e-12));
    }
}

TEST_CASE("cyclic_thomas_solve: 5x5 cyclic system", "[trajectory][tridiagonal]")
{
    constexpr int n = 5;
    std::vector<double> a(n), b(n), c(n), d(n);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd rhs(n);

    for (int i = 0; i < n; ++i) {
        b[static_cast<std::size_t>(i)] = 6.0;
        A(i, i) = 6.0;
        if (i > 0) {
            a[static_cast<std::size_t>(i)] = 1.0;
            A(i, i - 1) = 1.0;
        }
        if (i < n - 1) {
            c[static_cast<std::size_t>(i)] = 1.0;
            A(i, i + 1) = 1.0;
        }
        d[static_cast<std::size_t>(i)] = static_cast<double>(i * i + 1);
        rhs(i) = static_cast<double>(i * i + 1);
    }

    double alpha = 0.8;  // A[0][n-1]
    double beta = 0.3;   // A[n-1][0]
    A(0, n - 1) = alpha;
    A(n - 1, 0) = beta;

    Eigen::VectorXd x_ref = A.colPivHouseholderQr().solve(rhs);

    auto const solved = ctrlpp::detail::cyclic_thomas_solve(a, b, c, d, alpha, beta);
    REQUIRE(solved.has_value());

    for (int i = 0; i < n; ++i) {
        REQUIRE_THAT(d[static_cast<std::size_t>(i)], WithinAbs(x_ref(i), 1e-12));
    }
}

TEST_CASE("thomas_solve: a zero leading pivot is refused", "[trajectory][tridiagonal]")
{
    auto const solved = solve_tridiagonal({0.0, 1.0, 1.0}, {0.0, 2.0, 2.0}, {1.0, 1.0, 0.0});
    REQUIRE_FALSE(solved.has_value());
    REQUIRE(solved.error() == ctrlpp::spline_error::unrepresentable_spline);
}

TEST_CASE("thomas_solve: a zero 1x1 system is refused", "[trajectory][tridiagonal]")
{
    auto const solved = solve_tridiagonal({0.0}, {0.0}, {0.0});
    REQUIRE_FALSE(solved.has_value());
    REQUIRE(solved.error() == ctrlpp::spline_error::unrepresentable_spline);
}

TEST_CASE("thomas_solve: a NaN pivot is refused", "[trajectory][tridiagonal]")
{
    auto const nan = std::numeric_limits<double>::quiet_NaN();
    auto const solved = solve_tridiagonal({0.0, 1.0}, {nan, 1.0}, {1.0, 0.0});
    REQUIRE_FALSE(solved.has_value());
    REQUIRE(solved.error() == ctrlpp::spline_error::unrepresentable_spline);
}

TEST_CASE("thomas_solve: an exactly singular 2x2 is refused", "[trajectory][tridiagonal]")
{
    auto const solved = solve_tridiagonal({0.0, 1.0}, {1.0, 1.0}, {1.0, 0.0});
    REQUIRE_FALSE(solved.has_value());
    REQUIRE(solved.error() == ctrlpp::spline_error::unrepresentable_spline);
}

TEST_CASE("thomas_solve: a reduced pivot is judged against the rounding bound",
          "[trajectory][tridiagonal]")
{
    auto const eps = std::numeric_limits<double>::epsilon();

    auto const within = solve_tridiagonal({0.0, 1.0}, {1.0, 1.0 + eps}, {1.0, 0.0});
    REQUIRE_FALSE(within.has_value());
    REQUIRE(within.error() == ctrlpp::spline_error::unrepresentable_spline);

    std::vector<double> b = {1.0, 1.0 + 8.0 * eps};
    std::vector<double> d = {1.0, 1.0};
    auto const beyond = ctrlpp::detail::thomas_solve({0.0, 1.0}, b, {1.0, 0.0}, d);
    REQUIRE(beyond.has_value());
    REQUIRE_THAT(d[0], WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(d[1], WithinAbs(0.0, 1e-12));
}

TEST_CASE("cyclic_thomas_solve: a singular circulant is refused", "[trajectory][tridiagonal]")
{
    std::vector<double> a = {0.0, -1.0, -1.0, -1.0};
    std::vector<double> b = {2.0, 2.0, 2.0, 2.0};
    std::vector<double> c = {-1.0, -1.0, -1.0, 0.0};
    std::vector<double> d = {1.0, 0.0, -1.0, 0.0};

    auto const solved = ctrlpp::detail::cyclic_thomas_solve(a, b, c, d, -1.0, -1.0);
    REQUIRE_FALSE(solved.has_value());
    REQUIRE(solved.error() == ctrlpp::spline_error::unrepresentable_spline);
}

// -- Cubic spline tests -------------------------------------------------------

// -- Natural BC ---------------------------------------------------------------

TEST_CASE("Natural BC: passes through all waypoints", "[trajectory][cubic_spline][natural]")
{
    auto spline = make_spline({
        .times = {0.0, 1.0, 2.0, 3.0, 4.0},
        .positions = {0.0, 1.0, 0.0, 1.0, 0.0},
        .bc = ctrlpp::boundary_condition::natural,
    });

    std::vector<double> ts = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> qs = {0.0, 1.0, 0.0, 1.0, 0.0};

    for (std::size_t i = 0; i < ts.size(); ++i) {
        auto const pt = spline.evaluate(ts[i]);
        REQUIRE_THAT(pt.position[0], WithinAbs(qs[i], 1e-12));
    }
}

TEST_CASE("Natural BC: zero acceleration at endpoints", "[trajectory][cubic_spline][natural]")
{
    auto spline = make_spline({
        .times = {0.0, 1.0, 2.0, 3.0, 4.0},
        .positions = {0.0, 1.0, 0.0, 1.0, 0.0},
        .bc = ctrlpp::boundary_condition::natural,
    });

    auto const p0 = spline.evaluate(0.0);
    auto const pT = spline.evaluate(4.0);

    REQUIRE_THAT(p0.acceleration[0], WithinAbs(0.0, 1e-10));
    REQUIRE_THAT(pT.acceleration[0], WithinAbs(0.0, 1e-10));
}

TEST_CASE("Natural BC: C2 continuity at interior knots", "[trajectory][cubic_spline][natural]")
{
    auto spline = make_spline({
        .times = {0.0, 1.0, 2.0, 3.0, 4.0},
        .positions = {0.0, 1.0, 0.0, 1.0, 0.0},
        .bc = ctrlpp::boundary_condition::natural,
    });

    double const eps = 1e-8;
    for (double t : {1.0, 2.0, 3.0}) {
        auto const left = spline.evaluate(t - eps);
        auto const right = spline.evaluate(t + eps);
        auto const at = spline.evaluate(t);

        REQUIRE_THAT(left.velocity[0], WithinAbs(at.velocity[0], 1e-5));
        REQUIRE_THAT(right.velocity[0], WithinAbs(at.velocity[0], 1e-5));

        REQUIRE_THAT(left.acceleration[0], WithinAbs(at.acceleration[0], 1e-4));
        REQUIRE_THAT(right.acceleration[0], WithinAbs(at.acceleration[0], 1e-4));
    }
}

TEST_CASE("Natural BC: duration returns t_n - t_0", "[trajectory][cubic_spline][natural]")
{
    auto spline = make_spline({
        .times = {0.0, 1.0, 2.0, 3.0, 4.0},
        .positions = {0.0, 1.0, 0.0, 1.0, 0.0},
        .bc = ctrlpp::boundary_condition::natural,
    });

    REQUIRE_THAT(spline.duration(), WithinAbs(4.0, 1e-12));
}

// -- Clamped BC ---------------------------------------------------------------

TEST_CASE("Clamped BC: respects endpoint velocities", "[trajectory][cubic_spline][clamped]")
{
    auto spline = make_spline({
        .times = {0.0, 1.0, 2.0},
        .positions = {0.0, 1.0, 0.5},
        .bc = ctrlpp::boundary_condition::clamped,
        .v0 = 1.0,
        .vn = -1.0,
    });

    auto const p0 = spline.evaluate(0.0);
    auto const pT = spline.evaluate(2.0);

    REQUIRE_THAT(p0.velocity[0], WithinAbs(1.0, 1e-10));
    REQUIRE_THAT(pT.velocity[0], WithinAbs(-1.0, 1e-10));
}

TEST_CASE("Clamped BC: passes through all waypoints", "[trajectory][cubic_spline][clamped]")
{
    auto spline = make_spline({
        .times = {0.0, 1.0, 2.0},
        .positions = {0.0, 1.0, 0.5},
        .bc = ctrlpp::boundary_condition::clamped,
        .v0 = 1.0,
        .vn = -1.0,
    });

    REQUIRE_THAT(spline.evaluate(0.0).position[0], WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(spline.evaluate(1.0).position[0], WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(spline.evaluate(2.0).position[0], WithinAbs(0.5, 1e-12));
}

// -- Periodic BC --------------------------------------------------------------

TEST_CASE("Periodic BC: matching velocity and acceleration at endpoints", "[trajectory][cubic_spline][periodic]")
{
    auto spline = make_spline({
        .times = {0.0, 1.0, 2.0, 3.0},
        .positions = {1.0, 2.0, 0.5, 1.0},
        .bc = ctrlpp::boundary_condition::periodic,
    });

    auto const p0 = spline.evaluate(0.0);
    auto const pT = spline.evaluate(3.0);

    REQUIRE_THAT(p0.velocity[0], WithinAbs(pT.velocity[0], 1e-10));
    REQUIRE_THAT(p0.acceleration[0], WithinAbs(pT.acceleration[0], 1e-10));
}

TEST_CASE("Periodic BC: passes through all waypoints", "[trajectory][cubic_spline][periodic]")
{
    auto spline = make_spline({
        .times = {0.0, 1.0, 2.0, 3.0},
        .positions = {1.0, 2.0, 0.5, 1.0},
        .bc = ctrlpp::boundary_condition::periodic,
    });

    REQUIRE_THAT(spline.evaluate(0.0).position[0], WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(spline.evaluate(1.0).position[0], WithinAbs(2.0, 1e-12));
    REQUIRE_THAT(spline.evaluate(2.0).position[0], WithinAbs(0.5, 1e-12));
    REQUIRE_THAT(spline.evaluate(3.0).position[0], WithinAbs(1.0, 1e-12));
}

// -- Validation ---------------------------------------------------------------
//
// Before the create conversion these invalid configurations were guarded
// only by assert(), so release builds accepted them silently and computed on
// garbage. Each case pins the specific rejection enumerator.

TEST_CASE("create rejects fewer than 2 waypoints", "[trajectory][cubic_spline][validation]")
{
    auto const result = ctrlpp::cubic_spline<double>::create({
        .times = {0.0},
        .positions = {1.0},
    });

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::too_few_points);
}

TEST_CASE("create rejects times/positions length mismatch", "[trajectory][cubic_spline][validation]")
{
    auto const result = ctrlpp::cubic_spline<double>::create({
        .times = {0.0, 1.0, 2.0},
        .positions = {0.0, 1.0},
    });

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::size_mismatch);
}

TEST_CASE("create rejects non-increasing knot times", "[trajectory][cubic_spline][validation]")
{
    auto const result = ctrlpp::cubic_spline<double>::create({
        .times = {0.0, 2.0, 1.0},
        .positions = {0.0, 1.0, 0.5},
    });

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::non_increasing_times);
}

TEST_CASE("create rejects periodic BC with 2 waypoints", "[trajectory][cubic_spline][validation][periodic]")
{
    // A periodic spline needs at least 3 waypoints: with only 2 the cyclic
    // system for the interior velocities is empty and the closed curve is
    // degenerate. The pre-conversion constructor accepted this input silently
    // in release builds and indexed out of range inside the cyclic solver.
    auto const result = ctrlpp::cubic_spline<double>::create({
        .times = {0.0, 1.0},
        .positions = {1.0, 1.0},
        .bc = ctrlpp::boundary_condition::periodic,
    });

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::periodic_too_few_points);
}

TEST_CASE("create accepts periodic BC with 3 waypoints", "[trajectory][cubic_spline][validation][periodic]")
{
    auto const result = ctrlpp::cubic_spline<double>::create({
        .times = {0.0, 1.0, 2.0},
        .positions = {1.0, 2.0, 1.0},
        .bc = ctrlpp::boundary_condition::periodic,
    });

    REQUIRE(result.has_value());

    auto const p0 = result->evaluate(0.0);
    auto const pT = result->evaluate(2.0);

    REQUIRE_THAT(p0.position[0], WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(pT.position[0], WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(p0.velocity[0], WithinAbs(pT.velocity[0], 1e-10));
    REQUIRE_THAT(p0.acceleration[0], WithinAbs(pT.acceleration[0], 1e-10));
}

TEST_CASE("create rejects periodic BC with mismatched endpoints", "[trajectory][cubic_spline][validation][periodic]")
{
    auto const result = ctrlpp::cubic_spline<double>::create({
        .times = {0.0, 1.0, 2.0, 3.0},
        .positions = {1.0, 2.0, 0.5, 1.5},
        .bc = ctrlpp::boundary_condition::periodic,
    });

    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::periodic_endpoint_mismatch);
}

TEST_CASE("create accepts a valid natural configuration", "[trajectory][cubic_spline][validation]")
{
    auto const result = ctrlpp::cubic_spline<double>::create({
        .times = {0.0, 1.0, 2.0, 3.0, 4.0},
        .positions = {0.0, 1.0, 0.0, 1.0, 0.0},
        .bc = ctrlpp::boundary_condition::natural,
    });

    REQUIRE(result.has_value());

    std::vector<double> const ts = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> const qs = {0.0, 1.0, 0.0, 1.0, 0.0};
    for (std::size_t i = 0; i < ts.size(); ++i) {
        REQUIRE_THAT(result->evaluate(ts[i]).position[0], WithinAbs(qs[i], 1e-12));
    }
}

// -- Clamping behavior --------------------------------------------------------

TEST_CASE("All BCs: evaluate clamps outside [t0, tn]", "[trajectory][cubic_spline]")
{
    auto spline = make_spline({
        .times = {0.0, 1.0, 2.0},
        .positions = {0.0, 1.0, 0.5},
        .bc = ctrlpp::boundary_condition::natural,
    });

    auto const before = spline.evaluate(-1.0);
    auto const at0 = spline.evaluate(0.0);
    auto const after = spline.evaluate(10.0);
    auto const atT = spline.evaluate(2.0);

    REQUIRE_THAT(before.position[0], WithinAbs(at0.position[0], 1e-12));
    REQUIRE_THAT(after.position[0], WithinAbs(atT.position[0], 1e-12));
}
