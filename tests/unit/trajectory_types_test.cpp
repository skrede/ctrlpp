#include "ctrlpp/traj/trajectory_types.h"
#include "ctrlpp/traj/trajectory_segment.h"
#include "ctrlpp/traj/detail/polynomial_eval.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cstddef>

using namespace ctrlpp;
using Catch::Matchers::WithinAbs;

// --- trajectory_point tests ---

TEST_CASE("trajectory_point<double, 1> has correct field types", "[traj][types]")
{
    trajectory_point<double, 1> pt{};
    pt.position[0] = 1.0;
    pt.velocity[0] = 2.0;
    pt.acceleration[0] = 3.0;

    CHECK_THAT(pt.position[0], WithinAbs(1.0, 1e-15));
    CHECK_THAT(pt.velocity[0], WithinAbs(2.0, 1e-15));
    CHECK_THAT(pt.acceleration[0], WithinAbs(3.0, 1e-15));
}

TEST_CASE("trajectory_point<double, 3> has 3D vector fields", "[traj][types]")
{
    trajectory_point<double, 3> pt{};
    pt.position = Vector<double, 3>{1.0, 2.0, 3.0};
    pt.velocity = Vector<double, 3>{4.0, 5.0, 6.0};
    pt.acceleration = Vector<double, 3>{7.0, 8.0, 9.0};

    CHECK_THAT(pt.position[2], WithinAbs(3.0, 1e-15));
    CHECK_THAT(pt.velocity[1], WithinAbs(5.0, 1e-15));
    CHECK_THAT(pt.acceleration[0], WithinAbs(7.0, 1e-15));
}

// --- path_point tests ---

TEST_CASE("path_point default-initializes to zero", "[traj][types]")
{
    path_point<double> pt{};
    CHECK_THAT(pt.q, WithinAbs(0.0, 1e-15));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-15));
    CHECK_THAT(pt.ddq, WithinAbs(0.0, 1e-15));
    CHECK_THAT(pt.dddq, WithinAbs(0.0, 1e-15));
}

TEST_CASE("path_point stores assigned values", "[traj][types]")
{
    path_point<double> pt{.q = 0.5, .dq = 1.5, .ddq = -3.0, .dddq = 12.0};
    CHECK_THAT(pt.q, WithinAbs(0.5, 1e-15));
    CHECK_THAT(pt.dq, WithinAbs(1.5, 1e-15));
    CHECK_THAT(pt.ddq, WithinAbs(-3.0, 1e-15));
    CHECK_THAT(pt.dddq, WithinAbs(12.0, 1e-15));
}

// --- trajectory_segment concept tests ---

namespace {

struct mock_segment
{
    auto evaluate(double /*t*/) const -> trajectory_point<double, 2>
    {
        return {};
    }
    auto duration() const -> double { return 1.0; }
};

struct missing_duration_segment
{
    auto evaluate(double /*t*/) const -> trajectory_point<double, 2>
    {
        return {};
    }
};

struct missing_evaluate_segment
{
    auto duration() const -> double { return 1.0; }
};

} // namespace

TEST_CASE("trajectory_segment concept accepts valid type", "[traj][concept]")
{
    static_assert(trajectory_segment<mock_segment, double, 2>);
    SUCCEED();
}

TEST_CASE("trajectory_segment concept rejects missing duration", "[traj][concept]")
{
    static_assert(!trajectory_segment<missing_duration_segment, double, 2>);
    SUCCEED();
}

TEST_CASE("trajectory_segment concept rejects missing evaluate", "[traj][concept]")
{
    static_assert(!trajectory_segment<missing_evaluate_segment, double, 2>);
    SUCCEED();
}

// --- Horner evaluation tests ---

TEST_CASE("horner_eval with [1, 2, 3] at tau=0.5 gives 2.75", "[traj][horner]")
{
    // p(tau) = 1 + 2*tau + 3*tau^2 = 1 + 1 + 0.75 = 2.75
    using V = Vector<double, 1>;
    std::array<V, 3> c{V{1.0}, V{2.0}, V{3.0}};

    auto result = detail::horner_eval<double, 1, 3>(c, 0.5);
    CHECK_THAT(result[0], WithinAbs(2.75, 1e-14));
}

TEST_CASE("horner_deriv1 with [1, 2, 3] at tau=0.5 gives 5.0", "[traj][horner]")
{
    // p'(tau) = 2 + 6*tau = 2 + 3 = 5.0
    using V = Vector<double, 1>;
    std::array<V, 3> c{V{1.0}, V{2.0}, V{3.0}};

    auto result = detail::horner_deriv1<double, 1, 3>(c, 0.5);
    CHECK_THAT(result[0], WithinAbs(5.0, 1e-14));
}

TEST_CASE("horner_deriv2 with [1, 2, 3] at tau=0.5 gives 6.0", "[traj][horner]")
{
    // p''(tau) = 2*3 = 6.0 (constant)
    using V = Vector<double, 1>;
    std::array<V, 3> c{V{1.0}, V{2.0}, V{3.0}};

    auto result = detail::horner_deriv2<double, 1, 3>(c, 0.5);
    CHECK_THAT(result[0], WithinAbs(6.0, 1e-14));
}

TEST_CASE("horner_eval with 4 coefficients (cubic) at tau=0.5", "[traj][horner]")
{
    // p(tau) = 1 + 0*tau + 3*tau^2 + (-2)*tau^3
    // = 1 + 0 + 0.75 + (-0.25) = 1.5
    using V = Vector<double, 1>;
    std::array<V, 4> c{V{1.0}, V{0.0}, V{3.0}, V{-2.0}};

    auto result = detail::horner_eval<double, 1, 4>(c, 0.5);
    CHECK_THAT(result[0], WithinAbs(1.5, 1e-14));
}

TEST_CASE("horner_deriv3 with 4 coefficients gives constant 3rd derivative", "[traj][horner]")
{
    // p(tau) = c0 + c1*tau + c2*tau^2 + c3*tau^3
    // p'''(tau) = 6*c3 = 6*(-2) = -12
    using V = Vector<double, 1>;
    std::array<V, 4> c{V{1.0}, V{0.0}, V{3.0}, V{-2.0}};

    auto result = detail::horner_deriv3<double, 1, 4>(c, 0.5);
    CHECK_THAT(result[0], WithinAbs(-12.0, 1e-14));
}
