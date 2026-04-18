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

foreach(bench IN LISTS INTERNAL_BENCHMARKS)
    add_executable(${bench} ${bench}.cpp)
    target_link_libraries(${bench} PRIVATE ctrlpp::ctrlpp nanobench)
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
