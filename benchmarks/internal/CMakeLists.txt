# Internal benchmarks -- one executable per component family

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

# These link only ctrlpp + nanobench. The ctrlpp bench logic is throw-free, but
# each TU embeds nanobench's implementation (ANKERL_NANOBENCH_IMPLEMENT), and
# nanobench throws internally with no disable-exceptions mode (unlike Catch2's
# CATCH_CONFIG_DISABLE_EXCEPTIONS). So -fno-exceptions cannot be honored here and
# CTRLPP_BENCH_NOEXCEPT_FLAGS stays empty until nanobench is compiled out of the
# bench TUs; the throw-free dogfood is enforced on the unit-test tree instead.
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

# Profile benches for CARE performance analysis; require ct_optcon; built under -fno-omit-frame-pointer
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
