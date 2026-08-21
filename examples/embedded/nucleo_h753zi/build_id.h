#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_BUILD_ID_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_BUILD_ID_H

#include <cstddef>

namespace ctrlpp
{

constexpr std::size_t kBuildIdTextLength = 41;

// Renders the running image's NT_GNU_BUILD_ID descriptor as lowercase hex. A
// note that is absent or not of the expected shape renders as the fixed text
// "unavailable", so a missing fingerprint reads as a named absence and never as
// a plausible value a host check could match against.
void render_build_id(char (&text)[kBuildIdTextLength]) noexcept;

}

#endif
