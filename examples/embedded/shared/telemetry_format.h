#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_TELEMETRY_FORMAT_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_TELEMETRY_FORMAT_H

#include <bit>
#include <span>
#include <array>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <concepts>
#include <string_view>
#include <type_traits>

namespace ctrlpp {

template<class T>
concept binary_float = std::floating_point<T> && std::numeric_limits<T>::is_iec559 && (sizeof(T) == sizeof(std::uint32_t) || sizeof(T) == sizeof(std::uint64_t));

template<binary_float T>
using bit_pattern_t = std::conditional_t<sizeof(T) == sizeof(std::uint64_t), std::uint64_t, std::uint32_t>;

template<binary_float T>
inline constexpr std::size_t kBitsTextLength = 2 * sizeof(T) + 1;

// Renders the object representation rather than the value, so values that
// compare equal but differ in bits (zero and negative zero) render differently.
// A hexadecimal-float conversion would also be lossless, but whether this leg's
// trimmed printf supports one is unverified, and this needs no C library at all.
template<binary_float T>
constexpr std::string_view render_bits(T value, std::span<char, kBitsTextLength<T>> text) noexcept
{
    constexpr std::string_view digits = "0123456789ABCDEF";
    auto pattern                      = std::bit_cast<bit_pattern_t<T>>(value);
    for(std::size_t i = text.size() - 1; i-- > 0;)
    {
        text[i] = digits[static_cast<std::size_t>(pattern & 0xFU)];
        pattern >>= 4U;
    }
    text.back() = '\0';
    return {text.data(), text.size() - 1};
}

enum class stream_tag : std::uint8_t
{
    meta,
    canary,
    exact,
    spread,
    verdict
};

inline constexpr std::array<stream_tag, 5> kStreamTags{stream_tag::meta, stream_tag::canary, stream_tag::exact, stream_tag::spread, stream_tag::verdict};

constexpr std::string_view tag_text(stream_tag tag) noexcept
{
    switch(tag)
    {
        case stream_tag::meta:
            return "[meta] ";
        case stream_tag::canary:
            return "[canary] ";
        case stream_tag::exact:
            return "[exact] ";
        case stream_tag::spread:
            return "[spread] ";
        case stream_tag::verdict:
            return "[verdict] ";
    }
    return {};
}

// Exact columns are asserted bitwise identical across runs; spread columns
// carry a stated spread bound instead, because flash wait states and fetch
// state move a cycle count on identical inputs.
enum class claim_class : std::uint8_t
{
    exact,
    spread
};

struct claim_field
{
    std::string_view name;
    claim_class claim;
};

// Each class block's header lists its fields in this order, and its first field
// keys the rows: the loop index for the exact block, the region for the spread.
inline constexpr std::array<claim_field, 16> kClaimFields{{{"k", claim_class::exact},
                                                           {"t_s", claim_class::exact},
                                                           {"x0", claim_class::exact},
                                                           {"x1", claim_class::exact},
                                                           {"u", claim_class::exact},
                                                           {"region", claim_class::spread},
                                                           {"regime", claim_class::spread},
                                                           {"n", claim_class::spread},
                                                           {"n_basis", claim_class::spread},
                                                           {"block", claim_class::spread},
                                                           {"raw_min", claim_class::spread},
                                                           {"raw_median", claim_class::spread},
                                                           {"raw_max", claim_class::spread},
                                                           {"max_at", claim_class::spread},
                                                           {"per_iteration_raw_max", claim_class::spread},
                                                           {"wraps", claim_class::spread}}};

constexpr std::size_t column_count(claim_class claim) noexcept
{
    std::size_t count = 0;
    for(const claim_field &field : kClaimFields)
        count += field.claim == claim ? 1U : 0U;
    return count;
}

constexpr stream_tag tag_of(claim_class claim) noexcept
{
    return claim == claim_class::exact ? stream_tag::exact : stream_tag::spread;
}

}

#endif
