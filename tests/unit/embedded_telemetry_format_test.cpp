#include "../../examples/embedded/shared/telemetry_format.h"

#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <span>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <cstddef>
#include <charconv>
#include <algorithm>
#include <string_view>
#include <system_error>

namespace {

template<class T>
std::string rendered(T value)
{
    std::array<char, ctrlpp::kBitsTextLength<T>> text{};
    return std::string{ctrlpp::render_bits(value, std::span{text})};
}

template<class T>
T parsed(const std::string &text)
{
    ctrlpp::bit_pattern_t<T> pattern{};
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), pattern, 16);
    REQUIRE(error == std::errc{});
    REQUIRE(end == text.data() + text.size());
    return std::bit_cast<T>(pattern);
}

template<class T>
void require_round_trip(T value)
{
    const std::string text = rendered(value);
    REQUIRE(text.size() == 2 * sizeof(T));
    REQUIRE(std::bit_cast<ctrlpp::bit_pattern_t<T>>(parsed<T>(text)) == std::bit_cast<ctrlpp::bit_pattern_t<T>>(value));
}

// The control demo's final state, as the board computed it: x0 after the last of
// its 201 steps.
constexpr double kDecayedState = 0x1.3909951dc4784p-16;

}

TEST_CASE("A known value renders as its binary64 bit pattern", "[embedded][telemetry]")
{
    REQUIRE(rendered(1.0) == "3FF0000000000000");
    REQUIRE(rendered(-2.0) == "C000000000000000");
    REQUIRE(rendered(kDecayedState) == "3EF3909951DC4784");
    REQUIRE(rendered(1.0F) == "3F800000");
}

TEST_CASE("Values one unit in the last place apart render differently", "[embedded][telemetry]")
{
    const double next = std::nextafter(kDecayedState, 1.0);
    REQUIRE(next != kDecayedState);
    REQUIRE(rendered(next) != rendered(kDecayedState));
    REQUIRE(rendered(std::nextafter(1.0F, 2.0F)) != rendered(1.0F));
}

TEST_CASE("Zero and negative zero render differently", "[embedded][telemetry]")
{
    REQUIRE(0.0 == -0.0);
    REQUIRE(rendered(0.0) == "0000000000000000");
    REQUIRE(rendered(-0.0) == "8000000000000000");
    REQUIRE(rendered(0.0F) != rendered(-0.0F));
}

TEST_CASE("The rendering round-trips across magnitudes at a fixed width", "[embedded][telemetry]")
{
    constexpr std::array<double, 7> doubles{0x1p-1022, 0x1.fffffffffffffp+1023, 0.0, -0.0, kDecayedState, -0x1.2ef26e6dc21a2p+3, 0x0.0000000000001p-1022};
    for(const double value : doubles)
        require_round_trip(value);
    constexpr std::array<float, 5> floats{0x1p-126F, 0x1.fffffep+127F, 0.0F, -0.0F, 0x1.39099ap-16F};
    for(const float value : floats)
        require_round_trip(value);
    REQUIRE(std::numeric_limits<double>::min() == 0x1p-1022);
    REQUIRE(std::numeric_limits<double>::max() == 0x1.fffffffffffffp+1023);
}

TEST_CASE("Every claim field maps to exactly one class", "[embedded][telemetry]")
{
    std::size_t exact  = 0;
    std::size_t spread = 0;
    for(const ctrlpp::claim_field &field : ctrlpp::kClaimFields)
    {
        REQUIRE_FALSE(field.name.empty());
        std::size_t named = 0;
        for(const ctrlpp::claim_field &other : ctrlpp::kClaimFields)
            named += other.name == field.name ? 1U : 0U;
        REQUIRE(named == 1);
        exact += field.claim == ctrlpp::claim_class::exact ? 1U : 0U;
        spread += field.claim == ctrlpp::claim_class::spread ? 1U : 0U;
    }
    REQUIRE(exact == ctrlpp::column_count(ctrlpp::claim_class::exact));
    REQUIRE(spread == ctrlpp::column_count(ctrlpp::claim_class::spread));
    REQUIRE(exact + spread == ctrlpp::kClaimFields.size());
}

TEST_CASE("Each class header is keyed by its first field", "[embedded][telemetry]")
{
    const auto first_of = [](ctrlpp::claim_class claim)
    { return std::find_if(ctrlpp::kClaimFields.begin(), ctrlpp::kClaimFields.end(), [claim](const ctrlpp::claim_field &field) { return field.claim == claim; })->name; };
    REQUIRE(first_of(ctrlpp::claim_class::exact) == "k");
    REQUIRE(first_of(ctrlpp::claim_class::spread) == "region");
}

TEST_CASE("The two classes travel under distinct tags among five", "[embedded][telemetry]")
{
    REQUIRE(ctrlpp::tag_of(ctrlpp::claim_class::exact) == ctrlpp::stream_tag::exact);
    REQUIRE(ctrlpp::tag_of(ctrlpp::claim_class::spread) == ctrlpp::stream_tag::spread);
    for(const ctrlpp::stream_tag tag : ctrlpp::kStreamTags)
    {
        const std::string_view text = ctrlpp::tag_text(tag);
        REQUIRE(text.size() > 3);
        REQUIRE(text.front() == '[');
        REQUIRE(text.substr(text.size() - 2) == "] ");
        std::size_t named = 0;
        for(const ctrlpp::stream_tag other : ctrlpp::kStreamTags)
            named += ctrlpp::tag_text(other) == text ? 1U : 0U;
        REQUIRE(named == 1);
    }
}
