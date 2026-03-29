# Drake competitive benchmarks
# Drake uses Bazel as its build system and has no CMake FetchContent support.
# No AUR package exists. Requires manual installation with CMake find_package
# support (pre-built tar.gz available for Ubuntu/macOS only).
# See benchmarks/README.md for installation instructions.
find_package(drake REQUIRED)

add_executable(bench_lqr_vs_drake bench_lqr_vs_drake.cpp)
target_link_libraries(bench_lqr_vs_drake PRIVATE
    ctrlpp::ctrlpp
    nanobench
    drake::drake)
