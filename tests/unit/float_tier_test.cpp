// Float-Scalar runtime tier: exercises tolerances that were hardcoded as
// fixed absolutes calibrated for double's ~2.22e-16 machine epsilon. At
// float's ~1.19e-7 machine epsilon, ordinary rounding noise routinely
// exceeds such a fixed threshold, so a value that should read as
// "numerically zero at this Scalar's own precision" instead reads as
// "significant", sending the code down the wrong branch.
//
// place.h's conjugate-pair check and biquad.h's near-singular DC-gain guard
// now both scale their tolerances by the Scalar's own machine epsilon, so
// the two cases below are active regression tests that pass at float.

#include "ctrlpp/dsp/biquad.h"
#include "ctrlpp/control/place.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <complex>
#include <cstddef>

using namespace ctrlpp;
using Catch::Matchers::WithinAbs;

TEST_CASE("place() accepts a numerically-conjugate float pole pair the way it does at double", "[float][anchor]")
{
    // The pole pair below is the same closed-form conjugate pair computed two
    // independent (mathematically identical) ways -- sin(theta) and
    // sqrt(1-cos^2(theta)) -- so any nonzero gap between them is pure
    // floating-point rounding, scaled by the Scalar's own machine epsilon.
    // At double that gap is ~1e-16 (well under the hardcoded 1e-12, so
    // validate_conjugate_pairs accepts the pair); at float it is ~5e-7
    // (comfortably above 1e-12), so the identical construction is spuriously
    // rejected as "not a conjugate pair" and place() refuses the design with
    // place_error::poles_not_conjugate_symmetric.
    auto poles_for = []<typename Scalar>(Scalar theta, Scalar wn) -> std::array<std::complex<Scalar>, 2>
    {
        Scalar sigma = -wn * std::cos(theta);
        Scalar im_a = wn * std::sin(theta);
        Scalar im_b = wn * std::sqrt(Scalar{1} - std::cos(theta) * std::cos(theta));
        return {std::complex<Scalar>{sigma, im_a}, std::complex<Scalar>{sigma, -im_b}};
    };

    Eigen::Matrix<float, 2, 2> A;
    A << 0.0f, 1.0f, -1.0f, -0.5f;
    Eigen::Matrix<float, 2, 1> B;
    B << 0.0f, 1.0f;

    auto poles_f = poles_for(1.2f, 5.0f);
    auto K_f = place<float, 2, 1>(A, B, poles_f);

    Eigen::Matrix<double, 2, 2> A_ref = A.cast<double>();
    Eigen::Matrix<double, 2, 1> B_ref = B.cast<double>();
    auto poles_d = poles_for(1.2, 5.0);
    auto K_ref = place<double, 2, 1>(A_ref, B_ref, poles_d);

    REQUIRE(K_ref.has_value());
    REQUIRE(K_f.has_value());

    const float eps = std::numeric_limits<float>::epsilon();
    const float tol = 2.0f * eps * (1.0f + static_cast<float>(K_ref->norm()));
    REQUIRE_THAT((*K_f)(0, 0), WithinAbs(static_cast<float>((*K_ref)(0, 0)), tol));
    REQUIRE_THAT((*K_f)(0, 1), WithinAbs(static_cast<float>((*K_ref)(0, 1)), tol));
}

TEST_CASE("biquad steady-state reset() at float matches the double reference near a singular DC gain", "[float][anchor]")
{
    // Coefficients constructed so 1+a1+a2 (the DC-gain denominator) is
    // exactly "ten machine epsilons" of the instantiated Scalar -- as close
    // to zero as that Scalar can meaningfully resolve. At double, ten
    // epsilons (~2.2e-15) is safely under the hardcoded 1e-12 guard, so
    // reset() correctly takes the "treat as singular, zero the state" branch
    // and the very first sample out of process() equals the reset value
    // exactly (b0=1). At float, ten epsilons (~1.2e-6) exceeds 1e-12, so the
    // guard does not fire and reset() instead divides by this noise-floor
    // denominator, producing a wildly incorrect steady state.
    auto denom_probe = []<typename Scalar>(Scalar reset_value) -> Scalar
    {
        biquad_coeffs<Scalar> c;
        c.b0 = Scalar{1};
        c.b1 = Scalar{0};
        c.b2 = Scalar{0};
        c.a1 = Scalar{-1} + Scalar{10} * std::numeric_limits<Scalar>::epsilon();
        c.a2 = Scalar{0};
        biquad<Scalar> filt(c);
        filt.reset(reset_value);
        return filt.process(reset_value);
    };

    float y_f = denom_probe(2.0f);
    double y_d = denom_probe(2.0);

    REQUIRE_THAT(y_d, WithinAbs(2.0, 1e-9));

    const float eps = std::numeric_limits<float>::epsilon();
    const float tol = 2.0f * eps * (1.0f + std::abs(static_cast<float>(y_d)));
    REQUIRE_THAT(y_f, WithinAbs(static_cast<float>(y_d), tol));
}
