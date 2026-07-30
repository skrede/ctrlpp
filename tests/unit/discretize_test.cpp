#include "ctrlpp/model/discretize.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <complex>

TEST_CASE("zoh discretize double integrator")
{
    ctrlpp::continuous_state_space<double, 2, 1, 1> sys{};
    sys.A(0, 0) = 0.0;
    sys.A(0, 1) = 1.0;
    sys.A(1, 0) = 0.0;
    sys.A(1, 1) = 0.0;
    sys.B(0, 0) = 0.0;
    sys.B(1, 0) = 1.0;
    sys.C(0, 0) = 1.0;
    sys.C(0, 1) = 0.0;
    sys.D(0, 0) = 0.0;

    double dt = 0.1;
    auto dsys = ctrlpp::discretize(ctrlpp::zoh{}, sys, dt);

    CHECK_THAT(dsys.A(0, 0), Catch::Matchers::WithinAbs(1.0, 1e-10));
    CHECK_THAT(dsys.A(0, 1), Catch::Matchers::WithinAbs(0.1, 1e-10));
    CHECK_THAT(dsys.A(1, 0), Catch::Matchers::WithinAbs(0.0, 1e-10));
    CHECK_THAT(dsys.A(1, 1), Catch::Matchers::WithinAbs(1.0, 1e-10));

    CHECK_THAT(dsys.B(0, 0), Catch::Matchers::WithinAbs(0.005, 1e-10));
    CHECK_THAT(dsys.B(1, 0), Catch::Matchers::WithinAbs(0.1, 1e-10));
}

TEST_CASE("zoh discretize preserves C and D")
{
    ctrlpp::continuous_state_space<double, 2, 1, 1> sys{};
    sys.A(0, 0) = 0.0;
    sys.A(0, 1) = 1.0;
    sys.A(1, 0) = 0.0;
    sys.A(1, 1) = 0.0;
    sys.B(0, 0) = 0.0;
    sys.B(1, 0) = 1.0;
    sys.C(0, 0) = 1.0;
    sys.C(0, 1) = 0.0;
    sys.D(0, 0) = 0.0;

    auto dsys = ctrlpp::discretize(ctrlpp::zoh{}, sys, 0.1);

    CHECK_THAT(dsys.C(0, 0), Catch::Matchers::WithinAbs(1.0, 1e-12));
    CHECK_THAT(dsys.C(0, 1), Catch::Matchers::WithinAbs(0.0, 1e-12));
    CHECK_THAT(dsys.D(0, 0), Catch::Matchers::WithinAbs(0.0, 1e-12));
}

// Backward error for a fixed-size direct solve (or product of resolvent-scaled terms) grows
// linearly with problem dimension (Trefethen & Bau, Numerical Linear Algebra, Lecture 15); NX
// bounds the accumulated rounding across the discretization formula's matrix products/inverses.
// The additive 1.0 floors the scale so a comparison against an exactly-zero reference still
// admits an eps-sized tolerance rather than zero.
namespace
{
template <typename Scalar, std::size_t NX>
Scalar discretize_tolerance(Scalar reference_magnitude)
{
    return static_cast<Scalar>(NX) * std::numeric_limits<Scalar>::epsilon() * (Scalar{1} + reference_magnitude);
}
}

TEST_CASE("tustin discretize maps continuous poles through the bilinear transform")
{
    // Upper-triangular A: eigenvalues are the diagonal entries s1 = -2, s2 = -5, and the
    // discretized Ad from a triangular A stays triangular, so its diagonal entries are the
    // discrete poles produced by the bilinear map.
    ctrlpp::continuous_state_space<double, 2, 1, 1> sys{};
    sys.A(0, 0) = -2.0;
    sys.A(0, 1) = 1.0;
    sys.A(1, 0) = 0.0;
    sys.A(1, 1) = -5.0;
    sys.B(0, 0) = 0.0;
    sys.B(1, 0) = 1.0;
    sys.C(0, 0) = 1.0;
    sys.C(0, 1) = 0.0;
    sys.D(0, 0) = 0.0;

    const double dt = 0.05;
    const double s1 = sys.A(0, 0);
    const double s2 = sys.A(1, 1);

    auto dsys = ctrlpp::discretize(ctrlpp::tustin{}, sys, dt);

    const double z1 = (1.0 + s1 * dt / 2.0) / (1.0 - s1 * dt / 2.0);
    const double z2 = (1.0 + s2 * dt / 2.0) / (1.0 - s2 * dt / 2.0);

    CHECK_THAT(dsys.A(0, 0), Catch::Matchers::WithinAbs(z1, discretize_tolerance<double, 2>(std::abs(z1))));
    CHECK_THAT(dsys.A(1, 1), Catch::Matchers::WithinAbs(z2, discretize_tolerance<double, 2>(std::abs(z2))));
    CHECK_THAT(dsys.A(1, 0), Catch::Matchers::WithinAbs(0.0, discretize_tolerance<double, 2>(0.0)));
}

TEST_CASE("forward euler discretize maps continuous poles through the first-order expansion")
{
    ctrlpp::continuous_state_space<double, 2, 1, 1> sys{};
    sys.A(0, 0) = -2.0;
    sys.A(0, 1) = 1.0;
    sys.A(1, 0) = 0.0;
    sys.A(1, 1) = -5.0;
    sys.B(0, 0) = 0.0;
    sys.B(1, 0) = 1.0;
    sys.C(0, 0) = 1.0;
    sys.C(0, 1) = 0.0;
    sys.D(0, 0) = 0.0;

    const double dt = 0.05;
    const double s1 = sys.A(0, 0);
    const double s2 = sys.A(1, 1);

    auto dsys = ctrlpp::discretize(ctrlpp::forward_euler{}, sys, dt);

    const double z1 = 1.0 + s1 * dt;
    const double z2 = 1.0 + s2 * dt;

    CHECK_THAT(dsys.A(0, 0), Catch::Matchers::WithinAbs(z1, discretize_tolerance<double, 2>(std::abs(z1))));
    CHECK_THAT(dsys.A(1, 1), Catch::Matchers::WithinAbs(z2, discretize_tolerance<double, 2>(std::abs(z2))));
    CHECK_THAT(dsys.A(0, 1), Catch::Matchers::WithinAbs(dt, discretize_tolerance<double, 2>(dt)));
}

TEST_CASE("backward euler discretize maps continuous poles through the resolvent")
{
    ctrlpp::continuous_state_space<double, 2, 1, 1> sys{};
    sys.A(0, 0) = -2.0;
    sys.A(0, 1) = 1.0;
    sys.A(1, 0) = 0.0;
    sys.A(1, 1) = -5.0;
    sys.B(0, 0) = 0.0;
    sys.B(1, 0) = 1.0;
    sys.C(0, 0) = 1.0;
    sys.C(0, 1) = 0.0;
    sys.D(0, 0) = 0.0;

    const double dt = 0.05;
    const double s1 = sys.A(0, 0);
    const double s2 = sys.A(1, 1);

    auto dsys = ctrlpp::discretize(ctrlpp::backward_euler{}, sys, dt);

    const double z1 = 1.0 / (1.0 - s1 * dt);
    const double z2 = 1.0 / (1.0 - s2 * dt);

    CHECK_THAT(dsys.A(0, 0), Catch::Matchers::WithinAbs(z1, discretize_tolerance<double, 2>(std::abs(z1))));
    CHECK_THAT(dsys.A(1, 1), Catch::Matchers::WithinAbs(z2, discretize_tolerance<double, 2>(std::abs(z2))));
}

TEST_CASE("tustin discretize matches the closed-form biproper scalar system")
{
    ctrlpp::continuous_state_space<double, 1, 1, 1> sys{};
    const double s = -3.0;
    const double b = 2.0;
    const double c = 1.5;
    const double d = 0.5;
    sys.A(0, 0) = s;
    sys.B(0, 0) = b;
    sys.C(0, 0) = c;
    sys.D(0, 0) = d;

    const double dt = 0.02;
    auto dsys = ctrlpp::discretize(ctrlpp::tustin{}, sys, dt);

    const double forward_resolvent_inverse = 1.0 / (1.0 - s * dt / 2.0);
    const double Ad = forward_resolvent_inverse * (1.0 + s * dt / 2.0);
    const double Bd = forward_resolvent_inverse * b * dt;
    const double Cd = c * forward_resolvent_inverse;
    const double Dd = d + Cd * b * dt / 2.0;

    CHECK_THAT(dsys.A(0, 0), Catch::Matchers::WithinAbs(Ad, discretize_tolerance<double, 1>(std::abs(Ad))));
    CHECK_THAT(dsys.B(0, 0), Catch::Matchers::WithinAbs(Bd, discretize_tolerance<double, 1>(std::abs(Bd))));
    CHECK_THAT(dsys.C(0, 0), Catch::Matchers::WithinAbs(Cd, discretize_tolerance<double, 1>(std::abs(Cd))));
    CHECK_THAT(dsys.D(0, 0), Catch::Matchers::WithinAbs(Dd, discretize_tolerance<double, 1>(std::abs(Dd))));
}

TEST_CASE("forward euler discretize matches the closed-form scalar system")
{
    ctrlpp::continuous_state_space<double, 1, 1, 1> sys{};
    const double s = -3.0;
    const double b = 2.0;
    const double c = 1.5;
    const double d = 0.5;
    sys.A(0, 0) = s;
    sys.B(0, 0) = b;
    sys.C(0, 0) = c;
    sys.D(0, 0) = d;

    const double dt = 0.02;
    auto dsys = ctrlpp::discretize(ctrlpp::forward_euler{}, sys, dt);

    const double Ad = 1.0 + s * dt;
    const double Bd = b * dt;

    CHECK_THAT(dsys.A(0, 0), Catch::Matchers::WithinAbs(Ad, discretize_tolerance<double, 1>(std::abs(Ad))));
    CHECK_THAT(dsys.B(0, 0), Catch::Matchers::WithinAbs(Bd, discretize_tolerance<double, 1>(std::abs(Bd))));
    CHECK_THAT(dsys.C(0, 0), Catch::Matchers::WithinAbs(c, discretize_tolerance<double, 1>(std::abs(c))));
    CHECK_THAT(dsys.D(0, 0), Catch::Matchers::WithinAbs(d, discretize_tolerance<double, 1>(std::abs(d))));
}

TEST_CASE("backward euler discretize matches the closed-form biproper scalar system")
{
    ctrlpp::continuous_state_space<double, 1, 1, 1> sys{};
    const double s = -3.0;
    const double b = 2.0;
    const double c = 1.5;
    const double d = 0.5;
    sys.A(0, 0) = s;
    sys.B(0, 0) = b;
    sys.C(0, 0) = c;
    sys.D(0, 0) = d;

    const double dt = 0.02;
    auto dsys = ctrlpp::discretize(ctrlpp::backward_euler{}, sys, dt);

    const double backward_resolvent_inverse = 1.0 / (1.0 - s * dt);
    const double Ad = backward_resolvent_inverse;
    const double Bd = backward_resolvent_inverse * b * dt;
    const double Cd = c * backward_resolvent_inverse;
    const double Dd = d + Cd * b * dt;

    CHECK_THAT(dsys.A(0, 0), Catch::Matchers::WithinAbs(Ad, discretize_tolerance<double, 1>(std::abs(Ad))));
    CHECK_THAT(dsys.B(0, 0), Catch::Matchers::WithinAbs(Bd, discretize_tolerance<double, 1>(std::abs(Bd))));
    CHECK_THAT(dsys.C(0, 0), Catch::Matchers::WithinAbs(Cd, discretize_tolerance<double, 1>(std::abs(Cd))));
    CHECK_THAT(dsys.D(0, 0), Catch::Matchers::WithinAbs(Dd, discretize_tolerance<double, 1>(std::abs(Dd))));
}

TEST_CASE("tustin prewarp discretize agrees with the plain bilinear map at the warped period")
{
    // Rescaling dt by dt_warp = (2/w_c) * tan(w_c*dt/2) before the bilinear map is the
    // definition of prewarping; discretize(tustin_prewarp{w_c}, sys, dt) must therefore equal
    // discretize(tustin{}, sys, dt_warp) exactly (both compute the same formula).
    ctrlpp::continuous_state_space<double, 1, 1, 1> sys{};
    sys.A(0, 0) = -4.0;
    sys.B(0, 0) = 1.0;
    sys.C(0, 0) = 1.0;
    sys.D(0, 0) = 0.0;

    const double dt = 0.02;
    const double w_c = 20.0;

    auto dsys_warp = ctrlpp::discretize(ctrlpp::tustin_prewarp<double>{w_c}, sys, dt);

    const double dt_warp = (2.0 / w_c) * std::tan(w_c * dt / 2.0);
    auto dsys_ref = ctrlpp::discretize(ctrlpp::tustin{}, sys, dt_warp);

    CHECK_THAT(dsys_warp.A(0, 0), Catch::Matchers::WithinAbs(dsys_ref.A(0, 0), discretize_tolerance<double, 1>(std::abs(dsys_ref.A(0, 0)))));
    CHECK_THAT(dsys_warp.B(0, 0), Catch::Matchers::WithinAbs(dsys_ref.B(0, 0), discretize_tolerance<double, 1>(std::abs(dsys_ref.B(0, 0)))));
}

TEST_CASE("tustin prewarp maps the critical frequency onto the continuous imaginary axis exactly")
{
    // The analytic reason prewarping works: substituting z = e^(j*w_c*dt) into the inverse
    // bilinear map s = (2/dt_warp) * (z-1)/(z+1) recovers s = j*w_c exactly, for any dt, when
    // dt_warp = (2/w_c) * tan(w_c*dt/2). This is the in-process analytic oracle for the
    // tustin_prewarp tag, independent of any state-space system.
    const double w_c = 12.0;
    const double dt = 0.03;
    const double dt_warp = (2.0 / w_c) * std::tan(w_c * dt / 2.0);

    const std::complex<double> z = std::exp(std::complex<double>(0.0, w_c * dt));
    const std::complex<double> s_mapped = (2.0 / dt_warp) * (z - 1.0) / (z + 1.0);

    const double tol = discretize_tolerance<double, 1>(w_c);
    CHECK_THAT(s_mapped.real(), Catch::Matchers::WithinAbs(0.0, tol));
    CHECK_THAT(s_mapped.imag(), Catch::Matchers::WithinAbs(w_c, tol));
}
