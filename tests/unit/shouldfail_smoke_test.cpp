// This proves that the vendored Catch2 v3.7.1 honors the [!shouldfail]
// special tag: a TEST_CASE tagged [!shouldfail] has its result inverted, so
// an assertion that fails is reported as an overall PASS. The test harness
// relies on this inversion to hold bug-exposing regression anchors green
// until the underlying defect is fixed and the tag is removed.

#include <catch2/catch_test_macros.hpp>

TEST_CASE("shouldfail tag inverts a failing assertion into a passing test case", "[harness][!shouldfail]")
{
    REQUIRE(1 == 2);
}
