# HPIPM comparison benchmarks

find_package(hpipm REQUIRED)
find_package(blasfeo REQUIRED)

include(FetchContent)
FetchContent_Declare(hpipm_cpp
    GIT_REPOSITORY https://github.com/mayataka/hpipm-cpp.git
    GIT_TAG main
    GIT_SHALLOW TRUE
    EXCLUDE_FROM_ALL)
FetchContent_MakeAvailable(hpipm_cpp)

add_executable(bench_qp_vs_hpipm bench_qp_vs_hpipm.cpp)

target_link_libraries(bench_qp_vs_hpipm PRIVATE
    ctrlpp::ctrlpp
    ctrlpp::osqp
    nanobench
    hpipm-cpp::hpipm-cpp
    hpipm
    blasfeo
)
