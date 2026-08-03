#include "ctrlpp/detail/riccati_solution.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <limits>

namespace
{

template <int N>
using square = Eigen::Matrix<double, N, N>;

// A Frobenius magnitude of an N-by-N matrix of a single repeated entry is
// N * |entry|, which is what the two end cases below compare against. Both are
// stated as that product rather than as a decimal so the case says what it
// expects rather than what it observed.
template <int N>
auto uniform_magnitude(double entry) -> double
{
    return static_cast<double>(N) * std::abs(entry);
}

}

TEST_CASE("A Riccati magnitude survives a sum of squares that leaves the top "
          "of the range",
          "[riccati][magnitude]")
{
    // The plain member forms the sum of squares directly, so it reports an
    // infinity for an operand every entry of which is finite. An acceptance
    // comparison holding two of those is not evidence.
    const double entry = 1.0e200;
    const square<2> operand = square<2>::Constant(entry);
    const double expected = uniform_magnitude<2>(entry);

    CHECK_FALSE(std::isfinite(operand.norm()));

    const auto resolved = ctrlpp::detail::resolve_magnitude(operand);
    REQUIRE(resolved.resolved);
    CHECK(std::isfinite(resolved.value));
    CHECK(std::abs(resolved.value - expected)
          <= std::numeric_limits<double>::epsilon() * expected);
}

TEST_CASE("A Riccati magnitude survives a sum of squares that leaves the "
          "bottom of the range",
          "[riccati][magnitude]")
{
    // The end the continuous solver actually fails at. The plain member reports
    // exactly zero for an operand carrying a perfectly representable answer, and
    // a residual scale of zero turns its acceptance test into `0 <= 0`.
    const double entry = 1.0e-200;
    const square<2> operand = square<2>::Constant(entry);
    const double expected = uniform_magnitude<2>(entry);

    CHECK(operand.norm() == 0.0);

    const auto resolved = ctrlpp::detail::resolve_magnitude(operand);
    REQUIRE(resolved.resolved);
    CHECK(resolved.value != 0.0);
    CHECK(std::abs(resolved.value - expected)
          <= std::numeric_limits<double>::epsilon() * expected);

    // The operand that separates the two surviving candidates: a single entry at
    // the smallest positive subnormal. A magnitude safe only against gradual
    // underflow of the whole sum still returns zero here.
    square<2> smallest = square<2>::Zero();
    smallest(1, 0) = std::numeric_limits<double>::denorm_min();

    CHECK(smallest.norm() == 0.0);

    const auto smallest_resolved =
        ctrlpp::detail::resolve_magnitude(smallest);
    REQUIRE(smallest_resolved.resolved);
    CHECK(smallest_resolved.value
          == std::numeric_limits<double>::denorm_min());
}

TEST_CASE("A Riccati magnitude reports the ordinary case unchanged",
          "[riccati][magnitude]")
{
    const square<4> operand = square<4>::Constant(1.0);
    const auto resolved = ctrlpp::detail::resolve_magnitude(operand);

    REQUIRE(resolved.resolved);
    CHECK(std::abs(resolved.value - operand.norm())
          <= std::numeric_limits<double>::epsilon() * operand.norm());

    const square<4> zero = square<4>::Zero();
    const auto resolved_zero = ctrlpp::detail::resolve_magnitude(zero);
    REQUIRE(resolved_zero.resolved);
    CHECK(resolved_zero.value == 0.0);
}

TEST_CASE("A Riccati magnitude is unresolved on a non-finite operand",
          "[riccati][magnitude]")
{
    square<2> infinite = square<2>::Constant(1.0);
    infinite(0, 1) = std::numeric_limits<double>::infinity();
    CHECK_FALSE(ctrlpp::detail::resolve_magnitude(infinite).resolved);
    CHECK_FALSE(ctrlpp::detail::resolve_largest_entry(infinite).resolved);

    square<2> not_a_number = square<2>::Constant(1.0);
    not_a_number(1, 1) = std::numeric_limits<double>::quiet_NaN();
    CHECK_FALSE(ctrlpp::detail::resolve_magnitude(not_a_number).resolved);
    CHECK_FALSE(ctrlpp::detail::resolve_largest_entry(not_a_number).resolved);

    // Finite entries whose magnitude is nonetheless not representable are
    // unresolved rather than infinite.
    const square<2> beyond =
        square<2>::Constant(std::numeric_limits<double>::max());
    CHECK_FALSE(ctrlpp::detail::resolve_magnitude(beyond).resolved);
    CHECK(ctrlpp::detail::resolve_largest_entry(beyond).resolved);
}

TEST_CASE("Combining Riccati magnitudes propagates unresolved",
          "[riccati][magnitude]")
{
    const ctrlpp::detail::resolved_magnitude<double> resolved{2.0, true};
    const ctrlpp::detail::resolved_magnitude<double> larger{5.0, true};
    const ctrlpp::detail::resolved_magnitude<double> unresolved{0.0, false};

    const auto combined =
        ctrlpp::detail::largest_magnitude({resolved, larger});
    REQUIRE(combined.resolved);
    CHECK(combined.value == 5.0);

    CHECK_FALSE(
        ctrlpp::detail::largest_magnitude({resolved, unresolved}).resolved);
    CHECK_FALSE(
        ctrlpp::detail::largest_magnitude({unresolved, larger}).resolved);

    CHECK_FALSE(ctrlpp::detail::scaled_magnitude(unresolved, 2.0).resolved);
    CHECK_FALSE(
        ctrlpp::detail::scaled_magnitude(
            resolved, std::numeric_limits<double>::infinity())
            .resolved);
    CHECK_FALSE(
        ctrlpp::detail::scaled_magnitude(
            {std::numeric_limits<double>::max(), true},
            std::numeric_limits<double>::max())
            .resolved);

    CHECK_FALSE(
        ctrlpp::detail::counted_rounding_bound(unresolved, 4).resolved);
    const auto bound = ctrlpp::detail::counted_rounding_bound(resolved, 4);
    REQUIRE(bound.resolved);
    CHECK(bound.value
          == 4.0 * std::numeric_limits<double>::epsilon() * resolved.value);
}

TEST_CASE("A Riccati magnitude comparison refuses an unresolved side",
          "[riccati][magnitude]")
{
    const ctrlpp::detail::resolved_magnitude<double> small{1.0, true};
    const ctrlpp::detail::resolved_magnitude<double> large{2.0, true};
    const ctrlpp::detail::resolved_magnitude<double> unresolved{0.0, false};

    CHECK(ctrlpp::detail::magnitude_within(small, large));
    CHECK(ctrlpp::detail::magnitude_within(small, small));
    CHECK_FALSE(ctrlpp::detail::magnitude_within(large, small));

    CHECK_FALSE(ctrlpp::detail::magnitude_within(unresolved, unresolved));
    CHECK_FALSE(ctrlpp::detail::magnitude_within(unresolved, large));
    CHECK_FALSE(ctrlpp::detail::magnitude_within(small, unresolved));

    // Two infinities do not compare equal-or-less through this helper. Reaching
    // it requires unresolved operands, which is exactly what the two matrices
    // below produce.
    const square<2> beyond =
        square<2>::Constant(std::numeric_limits<double>::max());
    CHECK_FALSE(ctrlpp::detail::magnitude_within(
        ctrlpp::detail::resolve_magnitude(beyond),
        ctrlpp::detail::resolve_magnitude(beyond)));
}
