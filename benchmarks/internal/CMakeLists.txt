set(INTERNAL_BENCHMARKS
    bench_pid
    bench_lqr
    bench_dare
    bench_kalman
    bench_ekf
    bench_ukf
    bench_mrac
    bench_l1
    bench_sysid
    bench_trajectory
    bench_dsp
    bench_so3
)

# Both branches are empty: each TU embeds nanobench's implementation, which throws internally and
# has no disable-exceptions mode.
if(NOT CTRLPP_TESTS_WITH_EXCEPTIONS AND NOT MSVC)
    set(CTRLPP_BENCH_NOEXCEPT_FLAGS "")
else()
    set(CTRLPP_BENCH_NOEXCEPT_FLAGS "")
endif()

foreach(bench IN LISTS INTERNAL_BENCHMARKS)
    add_executable(${bench} ${bench}.cpp)
    target_link_libraries(${bench} PRIVATE ctrlpp::ctrlpp nanobench)
    target_compile_options(${bench} PRIVATE ${CTRLPP_BENCH_NOEXCEPT_FLAGS})
endforeach()

# CARE profiling benches, built with frame pointers for the profiler's call stacks.
find_package(ct_optcon QUIET)
if(ct_optcon_FOUND)
    foreach(bench bench_care_vs_ct_nx16 bench_care_vs_ct_nx30 bench_care_methods)
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${bench}.cpp")
            add_executable(${bench} ${bench}.cpp)
            target_link_libraries(${bench} PRIVATE ctrlpp::ctrlpp nanobench ct_optcon)
            target_compile_options(${bench} PRIVATE -fno-omit-frame-pointer)
        endif()
    endforeach()
endif()
