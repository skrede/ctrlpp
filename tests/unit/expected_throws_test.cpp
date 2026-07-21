// The owned expected's value() accessor is checked: on an error state it calls
// on_bad_expected_access(), which throws bad_expected_access when exceptions are
// enabled and std::abort()s otherwise. CHECK_THROWS is unavailable under
// CATCH_CONFIG_DISABLE_EXCEPTIONS, so this assertion lives in the exceptions
// carve-out tree. Every throw-free part of the surface is covered by
// expected_test.cpp in the default -fno-exceptions tree.
#include "ctrlpp/expected.h"

#include <catch2/catch_test_macros.hpp>

#include <string>
#include <utility>

namespace
{

enum class fetch_error
{
    not_found,
    timed_out,
};

}

TEST_CASE("expected value throws on an error state", "[expected][throws]")
{
    ctrlpp::expected<int, fetch_error> bad{ctrlpp::unexpected(fetch_error::timed_out)};
    CHECK_THROWS(bad.value());
}

TEST_CASE("expected value throws on an rvalue error state", "[expected][throws]")
{
    ctrlpp::expected<std::string, fetch_error> bad{ctrlpp::unexpected(fetch_error::not_found)};
    CHECK_THROWS(std::move(bad).value());
}
