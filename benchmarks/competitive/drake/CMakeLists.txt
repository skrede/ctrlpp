# Drake competitive benchmarks
# Drake uses Bazel as its build system and has no CMake FetchContent support.
# No AUR package exists. Requires manual installation with CMake find_package
# support (pre-built tar.gz available for Ubuntu/macOS only).
find_package(drake REQUIRED)

# Plan 03 adds benchmark executables here
