#ifndef HPP_GUARD_CTRLPP_DETAIL_EXPECTED_ACCESS_H
#define HPP_GUARD_CTRLPP_DETAIL_EXPECTED_ACCESS_H

#include "ctrlpp/config.h"

#include <cstdlib>
#include <exception>

namespace ctrlpp::detail
{

struct unexpect_t
{
    explicit unexpect_t() = default;
};

inline constexpr unexpect_t unexpect{};

#if CTRLPP_HAS_EXCEPTIONS
struct bad_expected_access : std::exception
{
    const char* what() const noexcept override { return "bad ctrlpp::expected access"; }
};
#endif

[[noreturn]] inline void on_bad_expected_access()
{
#if CTRLPP_HAS_EXCEPTIONS
    throw bad_expected_access{};
#else
    std::abort();
#endif
}

}

#endif
