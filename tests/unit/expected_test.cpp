#include "ctrlpp/expected.h"

#include <catch2/catch_test_macros.hpp>

#include <string>
#include <utility>
#include <version>

#if defined(__cpp_lib_expected)
    #include <expected>
#endif

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

TEST_CASE("unexpected deduces its error type on both compilers", "[expected]")
{
    auto u = ctrlpp::unexpected(fetch_error::not_found);
    static_assert(std::is_same_v<decltype(u), ctrlpp::unexpected<fetch_error>>);
    CHECK(u.error() == fetch_error::not_found);
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

TEST_CASE("expected copy construction preserves both states", "[expected]")
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

TEST_CASE("expected move construction transfers a non-trivial payload", "[expected]")
{
    ctrlpp::expected<std::string, fetch_error> src{std::string("dare")};
    const ctrlpp::expected<std::string, fetch_error> dst = std::move(src);
    CHECK(*dst == "dare");
}

TEST_CASE("expected copy assignment crosses value and error states", "[expected]")
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

TEST_CASE("expected move assignment crosses value and error states", "[expected]")
{
    ctrlpp::expected<std::string, fetch_error> target{ctrlpp::unexpected(fetch_error::not_found)};

    target = ctrlpp::expected<std::string, fetch_error>{std::string("moved-in")};
    REQUIRE(target.has_value());
    CHECK(*target == "moved-in");

    target = ctrlpp::expected<std::string, fetch_error>{ctrlpp::unexpected(fetch_error::timed_out)};
    REQUIRE(!target.has_value());
    CHECK(target.error() == fetch_error::timed_out);
}

TEST_CASE("expected value_or returns the fallback on error", "[expected]")
{
    ctrlpp::expected<int, fetch_error> bad{ctrlpp::unexpected(fetch_error::not_found)};
    CHECK(bad.value_or(99) == 99);
    CHECK((ctrlpp::expected<int, fetch_error>{7}).value_or(99) == 7);
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

TEST_CASE("expected void specialization copies and assigns across states", "[expected]")
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

#if defined(__cpp_lib_expected)

TEST_CASE("expected converts from a std::expected value and error", "[expected][interop]")
{
    const std::expected<int, fetch_error> std_ok{42};
    const ctrlpp::expected<int, fetch_error> ok{std_ok};
    REQUIRE(ok.has_value());
    CHECK(*ok == 42);

    const std::expected<int, fetch_error> std_bad{std::unexpected(fetch_error::timed_out)};
    const ctrlpp::expected<int, fetch_error> bad{std_bad};
    REQUIRE(!bad.has_value());
    CHECK(bad.error() == fetch_error::timed_out);
}

TEST_CASE("expected converts to a std::expected value and error", "[expected][interop]")
{
    const ctrlpp::expected<int, fetch_error> ok{7};
    const auto std_ok = static_cast<std::expected<int, fetch_error>>(ok);
    REQUIRE(std_ok.has_value());
    CHECK(*std_ok == 7);

    const ctrlpp::expected<int, fetch_error> bad{ctrlpp::unexpected(fetch_error::not_found)};
    const auto std_bad = static_cast<std::expected<int, fetch_error>>(bad);
    REQUIRE(!std_bad.has_value());
    CHECK(std_bad.error() == fetch_error::not_found);
}

TEST_CASE("void expected converts to and from std::expected", "[expected][interop]")
{
    const std::expected<void, fetch_error> std_bad{std::unexpected(fetch_error::timed_out)};
    const ctrlpp::expected<void, fetch_error> bad{std_bad};
    REQUIRE(!bad.has_value());
    CHECK(bad.error() == fetch_error::timed_out);

    const ctrlpp::expected<void, fetch_error> ok{};
    const auto std_ok = static_cast<std::expected<void, fetch_error>>(ok);
    CHECK(std_ok.has_value());
}

TEST_CASE("unexpected converts to and from std::unexpected", "[expected][interop]")
{
    const std::unexpected<fetch_error> std_u{fetch_error::not_found};
    const ctrlpp::unexpected<fetch_error> u{std_u};
    CHECK(u.error() == fetch_error::not_found);

    const auto back = static_cast<std::unexpected<fetch_error>>(u);
    CHECK(back.error() == fetch_error::not_found);
}

#endif
