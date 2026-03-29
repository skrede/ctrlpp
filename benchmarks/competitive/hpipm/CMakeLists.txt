# HPIPM competitive benchmarks
# No FetchContent fallback -- compiled libraries, AUR available (yay -S hpipm blasfeo)
find_package(hpipm REQUIRED)
find_package(blasfeo REQUIRED)

# hpipm-cpp C++ wrapper (header-only, FetchContent from GitHub)
include(FetchContent)
FetchContent_Declare(hpipm_cpp
    GIT_REPOSITORY https://github.com/mayataka/hpipm-cpp.git
    GIT_SHALLOW    TRUE
    EXCLUDE_FROM_ALL)
FetchContent_MakeAvailable(hpipm_cpp)

# Plan 03 adds benchmark executables here
