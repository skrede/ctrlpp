#include "ctrlpp/expected.h"

#include <catch2/catch_test_macros.hpp>

#include <string>
#include <cstdint>
#include <utility>
#include <version>
#include <type_traits>

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

struct pinned
{
    pinned(std::int32_t identifier, double scale)
        : id{identifier}
        , gain{scale}
    {
    }

    pinned(const pinned&) = delete;
    pinned(pinned&&) = delete;

    std::int32_t id;
    double gain;
};

struct tallied
{
    explicit tallied(std::int32_t& destructions)
        : destroyed{destructions}
    {
    }

    tallied(const tallied&) = delete;
    tallied(tallied&&) = delete;

    ~tallied()
    {
        ++destroyed;
    }

    std::int32_t& destroyed;
};

struct guarded
{
    explicit guarded(std::int32_t identifier)
        : id{identifier}
    {
    }

    guarded(const guarded& other) noexcept(false)
        : id{other.id}
    {
    }

    guarded(guarded&& other) noexcept(false)
        : id{other.id}
    {
    }

    auto operator=(const guarded& other) noexcept(false) -> guarded&
    {
        id = other.id;
        return *this;
    }

    auto operator=(guarded&& other) noexcept(false) -> guarded&
    {
        id = other.id;
        return *this;
    }

    std::int32_t id;
};

auto make_pinned(std::int32_t id) -> ctrlpp::expected<pinned, fetch_error>
{
    if(id < 0)
        return ctrlpp::unexpected(fetch_error::not_found);
    return ctrlpp::expected<pinned, fetch_error>{std::in_place, id, 0.5};
}

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

TEST_CASE("expected constructs a value it can neither copy nor move in place", "[expected]")
{
    const auto e = make_pinned(7);

    REQUIRE(e.has_value());
    CHECK(e->id == 7);
    CHECK(e->gain == 0.5);

    const auto refused = make_pinned(-1);
    REQUIRE(!refused.has_value());
    CHECK(refused.error() == fetch_error::not_found);
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

// A value type whose copy and move may both throw sends every error-to-value assignment through the
// path that has to keep the old error recoverable, which is the one an exception-free build reduces.
TEST_CASE("expected assigns across states when the value type's copy and move may throw", "[expected]")
{
    static_assert(!std::is_nothrow_copy_constructible_v<guarded> && !std::is_nothrow_move_constructible_v<guarded>);
    using guarded_result = ctrlpp::expected<guarded, fetch_error>;

    guarded_result moved_into{ctrlpp::unexpected(fetch_error::not_found)};
    moved_into = guarded_result{guarded{1}};
    REQUIRE(moved_into.has_value());
    CHECK(moved_into->id == 1);

    const guarded_result source{guarded{2}};
    guarded_result copied_into{ctrlpp::unexpected(fetch_error::not_found)};
    copied_into = source;
    REQUIRE(copied_into.has_value());
    CHECK(copied_into->id == 2);

    const guarded value{3};
    guarded_result assigned{ctrlpp::unexpected(fetch_error::timed_out)};
    assigned = value;
    REQUIRE(assigned.has_value());
    CHECK(assigned->id == 3);

    assigned = ctrlpp::unexpected(fetch_error::not_found);
    REQUIRE(!assigned.has_value());
    CHECK(assigned.error() == fetch_error::not_found);
}

TEST_CASE("expected assigns an error over a value it can neither copy nor move", "[expected]")
{
    static_assert(!std::is_copy_constructible_v<pinned> && !std::is_move_constructible_v<pinned>);
    static_assert(!std::is_copy_assignable_v<pinned> && !std::is_move_assignable_v<pinned>);

    ctrlpp::expected<pinned, fetch_error> moved_in{std::in_place, 7, 0.5};
    moved_in = ctrlpp::unexpected(fetch_error::timed_out);
    REQUIRE(!moved_in.has_value());
    CHECK(moved_in.error() == fetch_error::timed_out);

    const auto refusal = ctrlpp::unexpected(fetch_error::not_found);
    ctrlpp::expected<pinned, fetch_error> copied_in{std::in_place, 8, 1.5};
    copied_in = refusal;
    REQUIRE(!copied_in.has_value());
    CHECK(copied_in.error() == fetch_error::not_found);
}

TEST_CASE("expected destroys the value exactly once when an error is assigned over it", "[expected]")
{
    std::int32_t destructions = 0;
    {
        ctrlpp::expected<tallied, fetch_error> result{std::in_place, destructions};
        result = ctrlpp::unexpected(fetch_error::not_found);
        CHECK(destructions == 1);
        REQUIRE(!result.has_value());

        result = ctrlpp::unexpected(fetch_error::timed_out);
        CHECK(destructions == 1);
        CHECK(result.error() == fetch_error::timed_out);
    }
    CHECK(destructions == 1);
}

TEST_CASE("expected replaces the error when an error is assigned over an error", "[expected]")
{
    ctrlpp::expected<pinned, fetch_error> result{ctrlpp::unexpect, fetch_error::not_found};

    result = ctrlpp::unexpected(fetch_error::timed_out);
    REQUIRE(!result.has_value());
    CHECK(result.error() == fetch_error::timed_out);

    const auto refusal = ctrlpp::unexpected(fetch_error::not_found);
    result = refusal;
    REQUIRE(!result.has_value());
    CHECK(result.error() == fetch_error::not_found);
}

TEST_CASE("expected assigns an unexpected whose error converts to its error type", "[expected]")
{
    using described = ctrlpp::expected<std::int32_t, std::string>;
    static_assert(!std::is_assignable_v<described&, ctrlpp::unexpected<payload>>);

    described result{3};
    result = ctrlpp::unexpected("setup refused");
    REQUIRE(!result.has_value());
    CHECK(result.error() == "setup refused");

    const ctrlpp::unexpected<const char*> again{"setup refused twice"};
    result = again;
    REQUIRE(!result.has_value());
    CHECK(result.error() == "setup refused twice");
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

TEST_CASE("expected void specialization assigns an error over success and over an error", "[expected]")
{
    ctrlpp::expected<void, fetch_error> result{};

    result = ctrlpp::unexpected(fetch_error::timed_out);
    REQUIRE(!result.has_value());
    CHECK(result.error() == fetch_error::timed_out);

    const auto refusal = ctrlpp::unexpected(fetch_error::not_found);
    result = refusal;
    REQUIRE(!result.has_value());
    CHECK(result.error() == fetch_error::not_found);

    ctrlpp::expected<void, fetch_error> copied_in{};
    copied_in = refusal;
    REQUIRE(!copied_in.has_value());
    CHECK(copied_in.error() == fetch_error::not_found);

    copied_in = ctrlpp::unexpected(fetch_error::timed_out);
    REQUIRE(!copied_in.has_value());
    CHECK(copied_in.error() == fetch_error::timed_out);

    ctrlpp::expected<void, std::string> described{};
    described = ctrlpp::unexpected("setup refused");
    REQUIRE(!described.has_value());
    CHECK(described.error() == "setup refused");

    const ctrlpp::unexpected<const char*> again{"setup refused twice"};
    described = again;
    REQUIRE(!described.has_value());
    CHECK(described.error() == "setup refused twice");
}

TEST_CASE("expected exposes the std member types on both specializations", "[expected]")
{
    using valued = ctrlpp::expected<pinned, fetch_error>;
    static_assert(std::is_same_v<valued::value_type, pinned>);
    static_assert(std::is_same_v<valued::error_type, fetch_error>);
    static_assert(std::is_same_v<valued::unexpected_type, ctrlpp::unexpected<fetch_error>>);

    using empty = ctrlpp::expected<void, std::string>;
    static_assert(std::is_same_v<empty::value_type, void>);
    static_assert(std::is_same_v<empty::error_type, std::string>);
    static_assert(std::is_same_v<empty::unexpected_type, ctrlpp::unexpected<std::string>>);
    SUCCEED();
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
