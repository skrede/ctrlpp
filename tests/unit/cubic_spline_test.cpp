#include "ctrlpp/traj/cubic_spline.h"
#include "ctrlpp/traj/detail/tridiagonal.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>
#include <cmath>
#include <vector>

using Catch::Matchers::WithinAbs;

// -- Tridiagonal solver tests -------------------------------------------------

TEST_CASE("thomas_solve: 3x3 system", "[traj][tridiagonal]")
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

    ctrlpp::detail::thomas_solve(a, b, c, d);

    for (int i = 0; i < 3; ++i) {
        REQUIRE_THAT(d[static_cast<std::size_t>(i)], WithinAbs(x_ref(i), 1e-12));
    }
}

TEST_CASE("thomas_solve: 5x5 diagonal-dominant system", "[traj][tridiagonal]")
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

    ctrlpp::detail::thomas_solve(a, b, c, d);

    for (int i = 0; i < n; ++i) {
        REQUIRE_THAT(d[static_cast<std::size_t>(i)], WithinAbs(x_ref(i), 1e-12));
    }
}

TEST_CASE("cyclic_thomas_solve: 3x3 cyclic system", "[traj][tridiagonal]")
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

    ctrlpp::detail::cyclic_thomas_solve(a, b, c, d, alpha, beta);

    for (int i = 0; i < n; ++i) {
        REQUIRE_THAT(d[static_cast<std::size_t>(i)], WithinAbs(x_ref(i), 1e-12));
    }
}

TEST_CASE("cyclic_thomas_solve: 5x5 cyclic system", "[traj][tridiagonal]")
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

    ctrlpp::detail::cyclic_thomas_solve(a, b, c, d, alpha, beta);

    for (int i = 0; i < n; ++i) {
        REQUIRE_THAT(d[static_cast<std::size_t>(i)], WithinAbs(x_ref(i), 1e-12));
    }
}

// -- Cubic spline tests -------------------------------------------------------

// -- Natural BC ---------------------------------------------------------------

TEST_CASE("Natural BC: passes through all waypoints", "[traj][cubic_spline][natural]")
{
    ctrlpp::cubic_spline<double> spline({
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

TEST_CASE("Natural BC: zero acceleration at endpoints", "[traj][cubic_spline][natural]")
{
    ctrlpp::cubic_spline<double> spline({
        .times = {0.0, 1.0, 2.0, 3.0, 4.0},
        .positions = {0.0, 1.0, 0.0, 1.0, 0.0},
        .bc = ctrlpp::boundary_condition::natural,
    });

    auto const p0 = spline.evaluate(0.0);
    auto const pT = spline.evaluate(4.0);

    REQUIRE_THAT(p0.acceleration[0], WithinAbs(0.0, 1e-10));
    REQUIRE_THAT(pT.acceleration[0], WithinAbs(0.0, 1e-10));
}

TEST_CASE("Natural BC: C2 continuity at interior knots", "[traj][cubic_spline][natural]")
{
    ctrlpp::cubic_spline<double> spline({
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

TEST_CASE("Natural BC: duration returns t_n - t_0", "[traj][cubic_spline][natural]")
{
    ctrlpp::cubic_spline<double> spline({
        .times = {0.0, 1.0, 2.0, 3.0, 4.0},
        .positions = {0.0, 1.0, 0.0, 1.0, 0.0},
        .bc = ctrlpp::boundary_condition::natural,
    });

    REQUIRE_THAT(spline.duration(), WithinAbs(4.0, 1e-12));
}

// -- Clamped BC ---------------------------------------------------------------

TEST_CASE("Clamped BC: respects endpoint velocities", "[traj][cubic_spline][clamped]")
{
    ctrlpp::cubic_spline<double> spline({
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

TEST_CASE("Clamped BC: passes through all waypoints", "[traj][cubic_spline][clamped]")
{
    ctrlpp::cubic_spline<double> spline({
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

TEST_CASE("Periodic BC: matching velocity and acceleration at endpoints", "[traj][cubic_spline][periodic]")
{
    ctrlpp::cubic_spline<double> spline({
        .times = {0.0, 1.0, 2.0, 3.0},
        .positions = {1.0, 2.0, 0.5, 1.0},
        .bc = ctrlpp::boundary_condition::periodic,
    });

    auto const p0 = spline.evaluate(0.0);
    auto const pT = spline.evaluate(3.0);

    REQUIRE_THAT(p0.velocity[0], WithinAbs(pT.velocity[0], 1e-10));
    REQUIRE_THAT(p0.acceleration[0], WithinAbs(pT.acceleration[0], 1e-10));
}

TEST_CASE("Periodic BC: passes through all waypoints", "[traj][cubic_spline][periodic]")
{
    ctrlpp::cubic_spline<double> spline({
        .times = {0.0, 1.0, 2.0, 3.0},
        .positions = {1.0, 2.0, 0.5, 1.0},
        .bc = ctrlpp::boundary_condition::periodic,
    });

    REQUIRE_THAT(spline.evaluate(0.0).position[0], WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(spline.evaluate(1.0).position[0], WithinAbs(2.0, 1e-12));
    REQUIRE_THAT(spline.evaluate(2.0).position[0], WithinAbs(0.5, 1e-12));
    REQUIRE_THAT(spline.evaluate(3.0).position[0], WithinAbs(1.0, 1e-12));
}

// -- Clamping behavior --------------------------------------------------------

TEST_CASE("All BCs: evaluate clamps outside [t0, tn]", "[traj][cubic_spline]")
{
    ctrlpp::cubic_spline<double> spline({
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
