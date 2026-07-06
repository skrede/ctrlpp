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

TEST_CASE("expected constructs from a value", "[expected]")
{
    ctrlpp::expected<int, fetch_error> e{42};

    REQUIRE(e.has_value());
    REQUIRE(static_cast<bool>(e));
    CHECK(*e == 42);
    CHECK(e.value() == 42);
}

TEST_CASE("expected constructs from a convertible value", "[expected]")
{
    ctrlpp::expected<std::string, fetch_error> e{"riccati"};

    REQUIRE(e.has_value());
    CHECK(*e == "riccati");
}

TEST_CASE("expected constructs from unexpected", "[expected]")
{
    ctrlpp::expected<int, fetch_error> e{ctrlpp::unexpected(fetch_error::not_found)};

    REQUIRE(!e.has_value());
    REQUIRE(!static_cast<bool>(e));
    CHECK(e.error() == fetch_error::not_found);
}

TEST_CASE("expected constructs the error in place via unexpect", "[expected]")
{
    ctrlpp::expected<int, fetch_error> e{ctrlpp::unexpect, fetch_error::timed_out};

    REQUIRE(!e.has_value());
    CHECK(e.error() == fetch_error::timed_out);
}

TEST_CASE("expected has_value and operator bool agree on both states", "[expected]")
{
    ctrlpp::expected<int, fetch_error> ok{7};
    ctrlpp::expected<int, fetch_error> bad{ctrlpp::unexpected(fetch_error::timed_out)};

    CHECK(ok.has_value());
    CHECK(static_cast<bool>(ok));
    CHECK(!bad.has_value());
    CHECK(!static_cast<bool>(bad));
}

TEST_CASE("expected exposes members through operator arrow", "[expected]")
{
    ctrlpp::expected<payload, fetch_error> e{payload{3, 0.5}};

    REQUIRE(e.has_value());
    CHECK(e->id == 3);
    CHECK(e->gain == 0.5);

    const ctrlpp::expected<payload, fetch_error> ce{payload{9, 1.5}};
    CHECK(ce->id == 9);
    CHECK((*ce).gain == 1.5);
}

TEST_CASE("expected value moves out of an rvalue", "[expected]")
{
    ctrlpp::expected<std::string, fetch_error> e{std::string("hamiltonian")};

    const std::string moved = std::move(e).value();
    CHECK(moved == "hamiltonian");
}

TEST_CASE("expected error is accessible on const and rvalue paths", "[expected]")
{
    const ctrlpp::expected<int, fetch_error> ce{ctrlpp::unexpected(fetch_error::not_found)};
    CHECK(ce.error() == fetch_error::not_found);

    ctrlpp::expected<int, fetch_error> e{ctrlpp::unexpected(fetch_error::timed_out)};
    CHECK(std::move(e).error() == fetch_error::timed_out);
}

TEST_CASE("expected void specialization reports success by default", "[expected]")
{
    ctrlpp::expected<void, fetch_error> e{};

    REQUIRE(e.has_value());
    REQUIRE(static_cast<bool>(e));
}

TEST_CASE("expected void specialization carries an error", "[expected]")
{
    ctrlpp::expected<void, fetch_error> e{ctrlpp::unexpected(fetch_error::timed_out)};

    REQUIRE(!e.has_value());
    REQUIRE(!static_cast<bool>(e));
    CHECK(e.error() == fetch_error::timed_out);
}

TEST_CASE("expected void specialization constructs the error in place", "[expected]")
{
    ctrlpp::expected<void, fetch_error> e{ctrlpp::unexpect, fetch_error::not_found};

    REQUIRE(!e.has_value());
    CHECK(e.error() == fetch_error::not_found);
}
