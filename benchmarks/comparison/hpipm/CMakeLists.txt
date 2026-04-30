# HPIPM comparison benchmark.
#
# Wires the HPIPM C API directly (the upstream hpipm-cpp wrapper does not
# build against current HPIPM struct definitions) and FetchContents both
# HPIPM and BLASFEO from upstream because the Arch hpipm-git / blasfeo-git
# packages ship LTO-only static archives with no usable symbol table.

include(FetchContent)

FetchContent_Declare(blasfeo
    GIT_REPOSITORY https://github.com/giaf/blasfeo.git
    GIT_TAG master
    GIT_SHALLOW TRUE
    EXCLUDE_FROM_ALL)

set(BLASFEO_TARGET GENERIC CACHE STRING "" FORCE)
set(LA HIGH_PERFORMANCE CACHE STRING "" FORCE)
set(BLASFEO_TESTING OFF CACHE BOOL "" FORCE)
set(BLASFEO_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BLASFEO_BENCHMARKS OFF CACHE BOOL "" FORCE)

# Both HPIPM 0.1.4.2 and BLASFEO 0.1.4.2 ship CMake 2.x minimums; lift to 3.5
# for this configure pass only.
set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "" FORCE)

FetchContent_MakeAvailable(blasfeo)

FetchContent_Declare(hpipm
    GIT_REPOSITORY https://github.com/giaf/hpipm.git
    GIT_TAG master
    GIT_SHALLOW TRUE
    EXCLUDE_FROM_ALL)

set(HPIPM_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(hpipm)

add_executable(bench_qp_vs_hpipm bench_qp_vs_hpipm.cpp)

target_include_directories(bench_qp_vs_hpipm PRIVATE
    ${blasfeo_SOURCE_DIR}/include
    ${hpipm_SOURCE_DIR}/include
)

target_link_libraries(bench_qp_vs_hpipm PRIVATE
    ctrlpp::ctrlpp
    ctrlpp::bench_problems
    ctrlpp::osqp
    nanobench
    hpipm
    blasfeo
)
