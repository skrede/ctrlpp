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
