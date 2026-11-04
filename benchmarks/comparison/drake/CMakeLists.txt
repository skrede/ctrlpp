# Drake builds with Bazel and has no FetchContent support; its pre-built archives cover Ubuntu and
# macOS only.
find_package(drake REQUIRED)

add_executable(bench_lqr_vs_drake bench_lqr_vs_drake.cpp)
target_link_libraries(bench_lqr_vs_drake PRIVATE
    ctrlpp::ctrlpp
    nanobench
    drake::drake
)
