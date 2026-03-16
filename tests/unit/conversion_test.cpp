#include "ctrlpp/conversion.h"
#include "ctrlpp/eigen_linalg.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using Policy = ctrlpp::EigenLinalgPolicy;

TEST_CASE("tf2ss converts first-order transfer function") {
    // H(s) = 1 / (s + 2)
    // Coefficients highest-degree-first: num = {1}, den = {1, 2}
    ctrlpp::TransferFunction<double, 0, 1, Policy> tf{
        .numerator = {1.0},
        .denominator = {1.0, 2.0}
    };

    auto ss = ctrlpp::tf2ss(tf);

    // Expected CCF: A = [-2], B = [1], C = [1], D = [0]
    CHECK_THAT(ss.A(0, 0), Catch::Matchers::WithinAbs(-2.0, 1e-12));
    CHECK_THAT(ss.B(0, 0), Catch::Matchers::WithinAbs(1.0, 1e-12));
    CHECK_THAT(ss.C(0, 0), Catch::Matchers::WithinAbs(1.0, 1e-12));
    CHECK_THAT(ss.D(0, 0), Catch::Matchers::WithinAbs(0.0, 1e-12));
}

TEST_CASE("tf2ss converts second-order transfer function") {
    // H(s) = (s + 1) / (s^2 + 3s + 2)
    // num = {1, 1}, den = {1, 3, 2}
    // NumDeg = 1, DenDeg = 2
    ctrlpp::TransferFunction<double, 1, 2, Policy> tf{
        .numerator = {1.0, 1.0},
        .denominator = {1.0, 3.0, 2.0}
    };

    auto ss = ctrlpp::tf2ss(tf);

    // CCF for 2nd order: A = [[0, 1], [-2, -3]], B = [[0], [1]], C = [1, 1], D = [0]
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
    // Known state-space in CCF for H(s) = 1/(s+2): A=[-2], B=[1], C=[1], D=[0]
    using SS = ctrlpp::ContinuousStateSpace<double, 1, 1, 1, Policy>;
    SS sys{};
    sys.A(0, 0) = -2.0;
    sys.B(0, 0) = 1.0;
    sys.C(0, 0) = 1.0;
    sys.D(0, 0) = 0.0;

    auto tf = ctrlpp::ss2tf(sys);

    // Expected: num = {0, 1} (0*s + 1), den = {1, 2} (s + 2)
    // Since ss2tf always returns NumDeg = NX = 1:
    // numer[0] = D = 0, numer[1] = C*I*B + D*p1 = 1 + 0*2 = 1
    CHECK_THAT(tf.numerator[0], Catch::Matchers::WithinAbs(0.0, 1e-12));
    CHECK_THAT(tf.numerator[1], Catch::Matchers::WithinAbs(1.0, 1e-12));
    CHECK_THAT(tf.denominator[0], Catch::Matchers::WithinAbs(1.0, 1e-12));
    CHECK_THAT(tf.denominator[1], Catch::Matchers::WithinAbs(2.0, 1e-12));
}

TEST_CASE("tf2ss -> ss2tf roundtrip preserves coefficients") {
    // H(s) = (s + 1) / (s^2 + 3s + 2)
    ctrlpp::TransferFunction<double, 1, 2, Policy> tf_original{
        .numerator = {1.0, 1.0},
        .denominator = {1.0, 3.0, 2.0}
    };

    auto ss = ctrlpp::tf2ss(tf_original);
    auto tf_recovered = ctrlpp::ss2tf(ss);

    // Denominator should match (monic)
    CHECK_THAT(tf_recovered.denominator[0], Catch::Matchers::WithinAbs(1.0, 1e-10));
    CHECK_THAT(tf_recovered.denominator[1], Catch::Matchers::WithinAbs(3.0, 1e-10));
    CHECK_THAT(tf_recovered.denominator[2], Catch::Matchers::WithinAbs(2.0, 1e-10));

    // Numerator: ss2tf returns NumDeg = NX = 2, so numer has 3 coefficients.
    // Original was degree 1 (= {1, 1}), so recovered should be {0, 1, 1}.
    CHECK_THAT(tf_recovered.numerator[0], Catch::Matchers::WithinAbs(0.0, 1e-10));
    CHECK_THAT(tf_recovered.numerator[1], Catch::Matchers::WithinAbs(1.0, 1e-10));
    CHECK_THAT(tf_recovered.numerator[2], Catch::Matchers::WithinAbs(1.0, 1e-10));
}

TEST_CASE("tf2ss handles proper transfer function with equal degrees") {
    // H(s) = (2s + 5) / (s + 3) -- NumDeg == DenDeg == 1
    ctrlpp::TransferFunction<double, 1, 1, Policy> tf{
        .numerator = {2.0, 5.0},
        .denominator = {1.0, 3.0}
    };

    auto ss = ctrlpp::tf2ss(tf);

    // CCF: A = [-3], B = [1], C = [5 - 2*3] = [-1], D = [2]
    CHECK_THAT(ss.A(0, 0), Catch::Matchers::WithinAbs(-3.0, 1e-12));
    CHECK_THAT(ss.B(0, 0), Catch::Matchers::WithinAbs(1.0, 1e-12));
    CHECK_THAT(ss.C(0, 0), Catch::Matchers::WithinAbs(-1.0, 1e-12));
    CHECK_THAT(ss.D(0, 0), Catch::Matchers::WithinAbs(2.0, 1e-12));
}
