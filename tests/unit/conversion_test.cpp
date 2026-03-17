#include "ctrlpp/conversion.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

TEST_CASE("tf2ss converts first-order transfer function") {
    ctrlpp::TransferFunction<double, 0, 1> tf{
        .numerator = {1.0},
        .denominator = {1.0, 2.0}
    };

    auto ss = ctrlpp::tf2ss(tf);

    CHECK_THAT(ss.A(0, 0), Catch::Matchers::WithinAbs(-2.0, 1e-12));
    CHECK_THAT(ss.B(0, 0), Catch::Matchers::WithinAbs(1.0, 1e-12));
    CHECK_THAT(ss.C(0, 0), Catch::Matchers::WithinAbs(1.0, 1e-12));
    CHECK_THAT(ss.D(0, 0), Catch::Matchers::WithinAbs(0.0, 1e-12));
}

TEST_CASE("tf2ss converts second-order transfer function") {
    ctrlpp::TransferFunction<double, 1, 2> tf{
        .numerator = {1.0, 1.0},
        .denominator = {1.0, 3.0, 2.0}
    };

    auto ss = ctrlpp::tf2ss(tf);

    CHECK_THAT(ss.A(0, 0), Catch::Matchers::WithinAbs(0.0, 1e-12));
    CHECK_THAT(ss.A(0, 1), Catch::Matchers::WithinAbs(1.0, 1e-12));
    CHECK_THAT(ss.A(1, 0), Catch::Matchers::WithinAbs(-2.0, 1e-12));
    CHECK_THAT(ss.A(1, 1), Catch::Matchers::WithinAbs(-3.0, 1e-12));

    CHECK_THAT(ss.B(0, 0), Catch::Matchers::WithinAbs(0.0, 1e-12));
    CHECK_THAT(ss.B(1, 0), Catch::Matchers::WithinAbs(1.0, 1e-12));

    CHECK_THAT(ss.C(0, 0), Catch::Matchers::WithinAbs(1.0, 1e-12));
    CHECK_THAT(ss.C(0, 1), Catch::Matchers::WithinAbs(1.0, 1e-12));

    CHECK_THAT(ss.D(0, 0), Catch::Matchers::WithinAbs(0.0, 1e-12));
}

TEST_CASE("ss2tf recovers transfer function coefficients") {
    using SS = ctrlpp::ContinuousStateSpace<double, 1, 1, 1>;
    SS sys{};
    sys.A(0, 0) = -2.0;
    sys.B(0, 0) = 1.0;
    sys.C(0, 0) = 1.0;
    sys.D(0, 0) = 0.0;

    auto tf = ctrlpp::ss2tf(sys);

    CHECK_THAT(tf.numerator[0], Catch::Matchers::WithinAbs(0.0, 1e-12));
    CHECK_THAT(tf.numerator[1], Catch::Matchers::WithinAbs(1.0, 1e-12));
    CHECK_THAT(tf.denominator[0], Catch::Matchers::WithinAbs(1.0, 1e-12));
    CHECK_THAT(tf.denominator[1], Catch::Matchers::WithinAbs(2.0, 1e-12));
}

TEST_CASE("tf2ss -> ss2tf roundtrip preserves coefficients") {
    ctrlpp::TransferFunction<double, 1, 2> tf_original{
        .numerator = {1.0, 1.0},
        .denominator = {1.0, 3.0, 2.0}
    };

    auto ss = ctrlpp::tf2ss(tf_original);
    auto tf_recovered = ctrlpp::ss2tf(ss);

    CHECK_THAT(tf_recovered.denominator[0], Catch::Matchers::WithinAbs(1.0, 1e-10));
    CHECK_THAT(tf_recovered.denominator[1], Catch::Matchers::WithinAbs(3.0, 1e-10));
    CHECK_THAT(tf_recovered.denominator[2], Catch::Matchers::WithinAbs(2.0, 1e-10));

    CHECK_THAT(tf_recovered.numerator[0], Catch::Matchers::WithinAbs(0.0, 1e-10));
    CHECK_THAT(tf_recovered.numerator[1], Catch::Matchers::WithinAbs(1.0, 1e-10));
    CHECK_THAT(tf_recovered.numerator[2], Catch::Matchers::WithinAbs(1.0, 1e-10));
}

TEST_CASE("tf2ss handles proper transfer function with equal degrees") {
    ctrlpp::TransferFunction<double, 1, 1> tf{
        .numerator = {2.0, 5.0},
        .denominator = {1.0, 3.0}
    };

    auto ss = ctrlpp::tf2ss(tf);

    CHECK_THAT(ss.A(0, 0), Catch::Matchers::WithinAbs(-3.0, 1e-12));
    CHECK_THAT(ss.B(0, 0), Catch::Matchers::WithinAbs(1.0, 1e-12));
    CHECK_THAT(ss.C(0, 0), Catch::Matchers::WithinAbs(-1.0, 1e-12));
    CHECK_THAT(ss.D(0, 0), Catch::Matchers::WithinAbs(2.0, 1e-12));
}
