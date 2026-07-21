// Exercises the hand-rolled C++20 fallback directly, even on a toolchain that
// ships <expected>. expected_test.cpp covers whichever target the switch
// selects (std::expected on modern toolchains); this TU pins the fallback so
// its raw-union storage, special members, and exception-gated value() stay
// under CI regardless of the host standard library.
#define CTRLPP_FORCE_EXPECTED_FALLBACK 1

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

struct payload
{
    int id;
    double gain;
};

}

TEST_CASE("fallback constructs from a value", "[expected][fallback]")
{
    ctrlpp::expected<int, fetch_error> e{42};

    REQUIRE(e.has_value());
    REQUIRE(static_cast<bool>(e));
    CHECK(*e == 42);
    CHECK(e.value() == 42);
}

TEST_CASE("fallback constructs from a convertible value", "[expected][fallback]")
{
    ctrlpp::expected<std::string, fetch_error> e{"riccati"};

    REQUIRE(e.has_value());
    CHECK(*e == "riccati");
}

TEST_CASE("fallback constructs from unexpected and in place", "[expected][fallback]")
{
    ctrlpp::expected<int, fetch_error> e{ctrlpp::unexpected(fetch_error::not_found)};
    REQUIRE(!e.has_value());
    CHECK(e.error() == fetch_error::not_found);

    ctrlpp::expected<int, fetch_error> in_place{ctrlpp::unexpect, fetch_error::timed_out};
    REQUIRE(!in_place.has_value());
    CHECK(in_place.error() == fetch_error::timed_out);
}

TEST_CASE("fallback exposes members through operator arrow and moves out", "[expected][fallback]")
{
    ctrlpp::expected<payload, fetch_error> e{payload{3, 0.5}};
    CHECK(e->id == 3);
    CHECK(e->gain == 0.5);

    ctrlpp::expected<std::string, fetch_error> s{std::string("hamiltonian")};
    const std::string moved = std::move(s).value();
    CHECK(moved == "hamiltonian");
}

TEST_CASE("fallback copy construction preserves both states", "[expected][fallback]")
{
    const ctrlpp::expected<std::string, fetch_error> ok{std::string("lqr")};
    const ctrlpp::expected<std::string, fetch_error> ok_copy = ok;
    CHECK(*ok_copy == "lqr");
    CHECK(*ok == "lqr");

    const ctrlpp::expected<std::string, fetch_error> bad{ctrlpp::unexpected(fetch_error::not_found)};
    const ctrlpp::expected<std::string, fetch_error> bad_copy = bad;
    CHECK(!bad_copy.has_value());
    CHECK(bad_copy.error() == fetch_error::not_found);
}

TEST_CASE("fallback move construction transfers a non-trivial payload", "[expected][fallback]")
{
    ctrlpp::expected<std::string, fetch_error> src{std::string("dare")};
    const ctrlpp::expected<std::string, fetch_error> dst = std::move(src);
    CHECK(*dst == "dare");
}

TEST_CASE("fallback copy assignment crosses value and error states", "[expected][fallback]")
{
    const ctrlpp::expected<std::string, fetch_error> value_src{std::string("kalman")};
    const ctrlpp::expected<std::string, fetch_error> error_src{ctrlpp::unexpected(fetch_error::timed_out)};

    ctrlpp::expected<std::string, fetch_error> target{std::string("seed")};

    target = error_src; // value -> error: destroy payload, construct error
    REQUIRE(!target.has_value());
    CHECK(target.error() == fetch_error::timed_out);

    target = value_src; // error -> value: destroy error, construct payload
    REQUIRE(target.has_value());
    CHECK(*target == "kalman");
}

TEST_CASE("fallback move assignment crosses value and error states", "[expected][fallback]")
{
    ctrlpp::expected<std::string, fetch_error> target{ctrlpp::unexpected(fetch_error::not_found)};

    target = ctrlpp::expected<std::string, fetch_error>{std::string("moved-in")};
    REQUIRE(target.has_value());
    CHECK(*target == "moved-in");

    target = ctrlpp::expected<std::string, fetch_error>{ctrlpp::unexpected(fetch_error::timed_out)};
    REQUIRE(!target.has_value());
    CHECK(target.error() == fetch_error::timed_out);
}

TEST_CASE("fallback value_or returns the fallback on error", "[expected][fallback]")
{
    ctrlpp::expected<int, fetch_error> bad{ctrlpp::unexpected(fetch_error::not_found)};
    CHECK(bad.value_or(99) == 99);
    CHECK((ctrlpp::expected<int, fetch_error>{7}).value_or(99) == 7);
}

TEST_CASE("fallback value throws on an error state", "[expected][fallback]")
{
    ctrlpp::expected<int, fetch_error> bad{ctrlpp::unexpected(fetch_error::timed_out)};
    CHECK_THROWS(bad.value());
}

TEST_CASE("fallback void specialization carries success and error", "[expected][fallback]")
{
    ctrlpp::expected<void, fetch_error> ok{};
    REQUIRE(ok.has_value());

    ctrlpp::expected<void, fetch_error> bad{ctrlpp::unexpected(fetch_error::timed_out)};
    REQUIRE(!bad.has_value());
    CHECK(bad.error() == fetch_error::timed_out);

    ctrlpp::expected<void, fetch_error> bad_in_place{ctrlpp::unexpect, fetch_error::not_found};
    CHECK(bad_in_place.error() == fetch_error::not_found);
}

TEST_CASE("fallback void specialization copies and assigns across states", "[expected][fallback]")
{
    const ctrlpp::expected<void, fetch_error> error_src{ctrlpp::unexpected(fetch_error::not_found)};
    const ctrlpp::expected<void, fetch_error> error_copy = error_src;
    CHECK(!error_copy.has_value());
    CHECK(error_copy.error() == fetch_error::not_found);

    ctrlpp::expected<void, fetch_error> target{};
    target = error_src; // value -> error
    REQUIRE(!target.has_value());
    CHECK(target.error() == fetch_error::not_found);

    target = ctrlpp::expected<void, fetch_error>{}; // error -> value
    CHECK(target.has_value());
}
