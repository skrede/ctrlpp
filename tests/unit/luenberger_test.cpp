#include "ctrlpp/estimation/luenberger.h"
#include "ctrlpp/control/place.h"

#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/model/state_space.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <complex>


TEST_CASE("luenberger observer convergence with known gain") {
    // Simple 2-state system: x(k+1) = A x(k) + B u(k), y(k) = C x(k)
    ctrlpp::discrete_state_space<double, 2, 1, 1> sys;
    sys.A << 0.9, 0.1,
             0.0, 0.8;
    sys.B << 0.0,
             1.0;
    sys.C << 1.0, 0.0;
    sys.D << 0.0;

    // Observer gain (manually chosen for stable A-LC)
    Eigen::Matrix<double, 2, 1> L;
    L << 0.5,
         0.3;

    // Initial estimate off from true state
    Eigen::Vector2d x0_est = Eigen::Vector2d::Zero();
    ctrlpp::luenberger_observer<double, 2, 1, 1> obs(sys, L, x0_est);

    // True state
    Eigen::Vector2d x_true;
    x_true << 1.0, 0.5;

    for (int i = 0; i < 30; ++i) {
        Eigen::Matrix<double, 1, 1> u;
        u << 0.0;

        // Simulate true system
        Eigen::Matrix<double, 1, 1> z = sys.C * x_true;
        x_true = sys.A * x_true + sys.B * u;

        // Observer predict-update
        obs.predict(u);
        obs.update(z);
    }

    auto est = obs.state();
    CHECK_THAT(est(0), Catch::Matchers::WithinAbs(x_true(0), 0.1));
    CHECK_THAT(est(1), Catch::Matchers::WithinAbs(x_true(1), 0.1));
}

TEST_CASE("place produces correct eigenvalues for 2-state SISO") {
    Eigen::Matrix<double, 2, 2> A;
    A << 0.0, 1.0,
        -2.0, -3.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.0,
         1.0;

    std::array<std::complex<double>, 2> desired = {
        std::complex<double>{-1.0, 1.0},
        std::complex<double>{-1.0, -1.0}
    };

    auto result = ctrlpp::place<double, 2, 1>(A, B, desired);
    REQUIRE(result.has_value());

    auto K = *result;
    Eigen::Matrix<double, 2, 2> Acl = A - B * K;

    // Check eigenvalues of A-BK match desired
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> solver(Acl, false);
    auto evals = solver.eigenvalues();

    // Sort by real then imaginary for comparison
    std::array<std::complex<double>, 2> computed = {evals(0), evals(1)};
    std::sort(computed.begin(), computed.end(),
        [](auto a, auto b) { return a.imag() < b.imag(); });
    std::sort(desired.begin(), desired.end(),
        [](auto a, auto b) { return a.imag() < b.imag(); });

    for (std::size_t i = 0; i < 2; ++i) {
        CHECK_THAT(computed[i].real(), Catch::Matchers::WithinAbs(desired[i].real(), 1e-8));
        CHECK_THAT(computed[i].imag(), Catch::Matchers::WithinAbs(desired[i].imag(), 1e-8));
    }
}

TEST_CASE("place_observer produces gain for convergent observer") {
    // Discrete system
    ctrlpp::discrete_state_space<double, 2, 1, 1> sys;
    sys.A << 1.0, 0.1,
             0.0, 1.0;
    sys.B << 0.005,
             0.1;
    sys.C << 1.0, 0.0;
    sys.D << 0.0;

    // Desired observer poles inside unit circle
    std::array<std::complex<double>, 2> desired_poles = {
        std::complex<double>{0.3, 0.0},
        std::complex<double>{0.2, 0.0}
    };

    auto L_opt = ctrlpp::place_observer<double, 2, 1>(sys.A, sys.C, desired_poles);
    REQUIRE(L_opt.has_value());

    auto L = *L_opt;
    Eigen::Vector2d x0_est = Eigen::Vector2d::Zero();
    ctrlpp::luenberger_observer<double, 2, 1, 1> obs(sys, L, x0_est);

    // True state
    Eigen::Vector2d x_true;
    x_true << 1.0, 0.5;

    for (int i = 0; i < 20; ++i) {
        Eigen::Matrix<double, 1, 1> u;
        u << 0.0;
        Eigen::Matrix<double, 1, 1> z = sys.C * x_true;
        x_true = sys.A * x_true + sys.B * u;

        obs.predict(u);
        obs.update(z);
    }

    auto est = obs.state();
    CHECK_THAT(est(0), Catch::Matchers::WithinAbs(x_true(0), 0.5));
    CHECK_THAT(est(1), Catch::Matchers::WithinAbs(x_true(1), 0.5));
}

TEST_CASE("place on uncontrollable system returns nullopt") {
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 0.0,
         0.0, 2.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.0,
         0.0;  // Zero B = uncontrollable

    std::array<std::complex<double>, 2> desired = {
        std::complex<double>{-1.0, 0.0},
        std::complex<double>{-2.0, 0.0}
    };

    auto result = ctrlpp::place<double, 2, 1>(A, B, desired);
    CHECK_FALSE(result.has_value());
}

TEST_CASE("luenberger MIMO observer with manual gain") {
    // NX=3, NU=1, NY=2
    ctrlpp::discrete_state_space<double, 3, 1, 2> sys;
    sys.A << 0.9, 0.1, 0.0,
             0.0, 0.8, 0.2,
             0.0, 0.0, 0.7;
    sys.B << 0.0,
             0.0,
             1.0;
    sys.C << 1.0, 0.0, 0.0,
             0.0, 1.0, 0.0;
    sys.D = Eigen::Matrix<double, 2, 1>::Zero();

    // Manual gain (3x2 matrix)
    Eigen::Matrix<double, 3, 2> L;
    L << 0.4, 0.0,
         0.0, 0.3,
         0.0, 0.1;

    Eigen::Vector3d x0_est = Eigen::Vector3d::Zero();
    ctrlpp::luenberger_observer<double, 3, 1, 2> obs(sys, L, x0_est);

    // Run a few cycles
    for (int i = 0; i < 10; ++i) {
        Eigen::Matrix<double, 1, 1> u;
        u << 0.1;
        obs.predict(u);

        Eigen::Vector2d z;
        z << 1.0, 0.5;
        obs.update(z);
    }

    auto est = obs.state();
    for (int i = 0; i < 3; ++i)
        CHECK(std::isfinite(est(i)));
}

TEST_CASE("luenberger set_model and set_gain") {
    ctrlpp::discrete_state_space<double, 2, 1, 1> sys;
    sys.A << 0.9, 0.1,
             0.0, 0.8;
    sys.B << 0.0,
             1.0;
    sys.C << 1.0, 0.0;
    sys.D << 0.0;

    Eigen::Matrix<double, 2, 1> L;
    L << 0.5, 0.3;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    ctrlpp::luenberger_observer<double, 2, 1, 1> obs(sys, L, x0);

    // Change gain
    Eigen::Matrix<double, 2, 1> L2;
    L2 << 0.8, 0.6;
    obs.set_gain(L2);

    // Change model
    auto sys2 = sys;
    sys2.A(0, 1) = 0.2;
    obs.set_model(sys2);

    // Reset state
    Eigen::Vector2d x_new;
    x_new << 1.0, 2.0;
    obs.reset(x_new);

    auto est = obs.state();
    CHECK_THAT(est(0), Catch::Matchers::WithinAbs(1.0, 1e-10));
    CHECK_THAT(est(1), Catch::Matchers::WithinAbs(2.0, 1e-10));
}

TEST_CASE("luenberger concept satisfaction") {
    static_assert(ctrlpp::ObserverPolicy<ctrlpp::luenberger_observer<double, 2, 1, 1>>);
    static_assert(!ctrlpp::CovarianceObserver<ctrlpp::luenberger_observer<double, 2, 1, 1>>);
    CHECK(true);
}

TEST_CASE("place with real poles") {
    Eigen::Matrix<double, 2, 2> A;
    A << 0.0, 1.0,
        -2.0, -3.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.0,
         1.0;

    std::array<std::complex<double>, 2> desired = {
        std::complex<double>{-5.0, 0.0},
        std::complex<double>{-6.0, 0.0}
    };

    auto result = ctrlpp::place<double, 2, 1>(A, B, desired);
    REQUIRE(result.has_value());

    auto K = *result;
    Eigen::Matrix<double, 2, 2> Acl = A - B * K;
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> solver(Acl, false);
    auto evals = solver.eigenvalues();

    std::array<double, 2> computed_real = {evals(0).real(), evals(1).real()};
    std::sort(computed_real.begin(), computed_real.end());
    CHECK_THAT(computed_real[0], Catch::Matchers::WithinAbs(-6.0, 1e-8));
    CHECK_THAT(computed_real[1], Catch::Matchers::WithinAbs(-5.0, 1e-8));
}
