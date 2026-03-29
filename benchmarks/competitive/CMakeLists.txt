option(CTRLPP_BENCH_LIBMPC "Build libmpc++ competitive benchmarks" OFF)
option(CTRLPP_BENCH_OSQP_EIGEN "Build osqp-eigen competitive benchmarks" OFF)
option(CTRLPP_BENCH_HPIPM "Build HPIPM competitive benchmarks" OFF)
option(CTRLPP_BENCH_CT "Build ETH Control Toolbox competitive benchmarks" OFF)
option(CTRLPP_BENCH_RUCKIG "Build ruckig competitive benchmarks" OFF)
option(CTRLPP_BENCH_DRAKE "Build Drake competitive benchmarks" OFF)

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
