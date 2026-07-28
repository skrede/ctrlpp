# Deliberately non-conforming fixture. Never configured by CMake, referenced by
# no add_subdirectory, part of no build.
#
# It exists so the rule asserting that the discard warning is promoted to an
# error can be shown to fail. The MSVC token is present and the GCC and Clang
# token is deliberately absent, so the rule reports exactly one violation here.

set(FIXTURE_WARNING_FLAGS
    $<$<CXX_COMPILER_ID:MSVC>:/W4 /we4834>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:
        -Wall
        -Wextra
    >
)
