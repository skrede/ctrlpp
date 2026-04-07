option(CTRLPP_BENCH_LIBMPC "Build libmpc++ comparison benchmarks" OFF)
option(CTRLPP_BENCH_OSQP_EIGEN "Build osqp-eigen comparison benchmarks" OFF)
option(CTRLPP_BENCH_HPIPM "Build HPIPM comparison benchmarks" OFF)
option(CTRLPP_BENCH_CT "Build ETH Control Toolbox comparison benchmarks" OFF)
option(CTRLPP_BENCH_RUCKIG "Build ruckig comparison benchmarks" OFF)
option(CTRLPP_BENCH_DRAKE "Build Drake comparison benchmarks" OFF)

if(CTRLPP_BENCH_LIBMPC)
    add_subdirectory(libmpc)
endif()

if(CTRLPP_BENCH_OSQP_EIGEN)
    add_subdirectory(osqp_eigen)
endif()

if(CTRLPP_BENCH_HPIPM)
    add_subdirectory(hpipm)
endif()

if(CTRLPP_BENCH_CT)
    add_subdirectory(ct)
endif()

if(CTRLPP_BENCH_RUCKIG)
    add_subdirectory(ruckig)
endif()

if(CTRLPP_BENCH_DRAKE)
    add_subdirectory(drake)
endif()

option(CTRLPP_BENCH_ARGMIN "Build nablapp comparison benchmarks" OFF)

if(CTRLPP_BENCH_ARGMIN)
    add_subdirectory(argmin)
endif()
