#include "ctrlpp/model/discretise.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEST_CASE("zoh discretise double integrator") {
    ctrlpp::continuous_state_space<double, 2, 1, 1> sys{};
    sys.A(0, 0) = 0.0; sys.A(0, 1) = 1.0;
    sys.A(1, 0) = 0.0; sys.A(1, 1) = 0.0;
    sys.B(0, 0) = 0.0; sys.B(1, 0) = 1.0;
    sys.C(0, 0) = 1.0; sys.C(0, 1) = 0.0;
    sys.D(0, 0) = 0.0;

    double dt = 0.1;
    auto dsys = ctrlpp::discretise(ctrlpp::zoh{}, sys, dt);

    CHECK_THAT(dsys.A(0, 0), Catch::Matchers::WithinAbs(1.0, 1e-10));
    CHECK_THAT(dsys.A(0, 1), Catch::Matchers::WithinAbs(0.1, 1e-10));
    CHECK_THAT(dsys.A(1, 0), Catch::Matchers::WithinAbs(0.0, 1e-10));
    CHECK_THAT(dsys.A(1, 1), Catch::Matchers::WithinAbs(1.0, 1e-10));

    CHECK_THAT(dsys.B(0, 0), Catch::Matchers::WithinAbs(0.005, 1e-10));
    CHECK_THAT(dsys.B(1, 0), Catch::Matchers::WithinAbs(0.1, 1e-10));
}

TEST_CASE("zoh discretise preserves C and D") {
    ctrlpp::continuous_state_space<double, 2, 1, 1> sys{};
    sys.A(0, 0) = 0.0; sys.A(0, 1) = 1.0;
    sys.A(1, 0) = 0.0; sys.A(1, 1) = 0.0;
    sys.B(0, 0) = 0.0; sys.B(1, 0) = 1.0;
    sys.C(0, 0) = 1.0; sys.C(0, 1) = 0.0;
    sys.D(0, 0) = 0.0;

    auto dsys = ctrlpp::discretise(ctrlpp::zoh{}, sys, 0.1);

    CHECK_THAT(dsys.C(0, 0), Catch::Matchers::WithinAbs(1.0, 1e-12));
    CHECK_THAT(dsys.C(0, 1), Catch::Matchers::WithinAbs(0.0, 1e-12));
    CHECK_THAT(dsys.D(0, 0), Catch::Matchers::WithinAbs(0.0, 1e-12));
}
