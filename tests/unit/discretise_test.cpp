#include "ctrlpp/eigen_discretise.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Policy = ctrlpp::EigenLinalgPolicy;

TEST_CASE("ZOH discretise double integrator") {
    // Continuous double integrator: x' = [[0,1],[0,0]]*x + [[0],[1]]*u
    ctrlpp::ContinuousStateSpace<Policy, double, 2, 1, 1> sys{};
    sys.A(0, 0) = 0.0; sys.A(0, 1) = 1.0;
    sys.A(1, 0) = 0.0; sys.A(1, 1) = 0.0;
    sys.B(0, 0) = 0.0; sys.B(1, 0) = 1.0;
    sys.C(0, 0) = 1.0; sys.C(0, 1) = 0.0;
    sys.D(0, 0) = 0.0;

    double dt = 0.1;
    auto dsys = ctrlpp::discretise(ctrlpp::ZOH{}, sys, dt);

    // Expected exact ZOH for double integrator:
    // Ad = [[1, dt], [0, 1]] = [[1, 0.1], [0, 1]]
    CHECK_THAT(dsys.A(0, 0), Catch::Matchers::WithinAbs(1.0, 1e-10));
    CHECK_THAT(dsys.A(0, 1), Catch::Matchers::WithinAbs(0.1, 1e-10));
    CHECK_THAT(dsys.A(1, 0), Catch::Matchers::WithinAbs(0.0, 1e-10));
    CHECK_THAT(dsys.A(1, 1), Catch::Matchers::WithinAbs(1.0, 1e-10));

    // Bd = [[dt^2/2], [dt]] = [[0.005], [0.1]]
    CHECK_THAT(dsys.B(0, 0), Catch::Matchers::WithinAbs(0.005, 1e-10));
    CHECK_THAT(dsys.B(1, 0), Catch::Matchers::WithinAbs(0.1, 1e-10));
}

TEST_CASE("ZOH discretise preserves C and D") {
    ctrlpp::ContinuousStateSpace<Policy, double, 2, 1, 1> sys{};
    sys.A(0, 0) = 0.0; sys.A(0, 1) = 1.0;
    sys.A(1, 0) = 0.0; sys.A(1, 1) = 0.0;
    sys.B(0, 0) = 0.0; sys.B(1, 0) = 1.0;
    sys.C(0, 0) = 1.0; sys.C(0, 1) = 0.0;
    sys.D(0, 0) = 0.0;

    auto dsys = ctrlpp::discretise(ctrlpp::ZOH{}, sys, 0.1);

    // Cd == C, Dd == D
    CHECK_THAT(dsys.C(0, 0), Catch::Matchers::WithinAbs(1.0, 1e-12));
    CHECK_THAT(dsys.C(0, 1), Catch::Matchers::WithinAbs(0.0, 1e-12));
    CHECK_THAT(dsys.D(0, 0), Catch::Matchers::WithinAbs(0.0, 1e-12));
}
