if(NOT CTRLPP_FETCH_BENCHMARK_DEPS)
    message(FATAL_ERROR "The cross-solver QP suite builds its solvers from pinned sources; re-run with -DCTRLPP_FETCH_BENCHMARK_DEPS=ON.")
endif()

include(${CMAKE_CURRENT_LIST_DIR}/../deps/argmin.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/../deps/qpoases.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/../deps/daqp.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/../deps/hpipm.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/../deps/tinympc.cmake)
find_package(proxsuite 0.7.3 EXACT REQUIRED)

# TinyMPC compiles against the Eigen it bundles (3.4.90), the rest of the suite against the system one,
# and both instantiate same-named Eigen templates whose definitions differ, which crashed the first
# build once the linker kept one copy for both. TinyMPC and the one unit that talks to it therefore
# live in a shared object that binds its own symbols and exports only the session interface.
add_library(ctrlpp_qp_solvers_tinympc SHARED tinympc_session.cpp)
target_link_libraries(ctrlpp_qp_solvers_tinympc PRIVATE tinympcstatic)
set_target_properties(ctrlpp_qp_solvers_tinympc PROPERTIES VISIBILITY_INLINES_HIDDEN ON)
target_link_options(ctrlpp_qp_solvers_tinympc PRIVATE LINKER:--exclude-libs,ALL LINKER:-Bsymbolic)

add_library(ctrlpp_qp_solvers STATIC
    argmin_arms.cpp
    argmin_settings.cpp
    osqp_arm.cpp
    qpoases_arm.cpp
    daqp_arm.cpp
    proxqp_arms.cpp
    hpipm_stages.cpp
    hpipm_arm.cpp
    hpipm_settings.cpp
    tinympc_arm.cpp
    arm_factory.cpp
    dual_recovery.cpp
    cases.cpp
    shaped_case.cpp
    reference.cpp
    binding_gate.cpp
    probe_problems.cpp
    probe.cpp
    instruction_counter.cpp
    frequency_log.cpp
    record_writer.cpp
    sequence.cpp
    blocks.cpp
    driver_options.cpp
    describe.cpp
    driver.cpp
)
target_include_directories(ctrlpp_qp_solvers PRIVATE
    ${qpoases_SOURCE_DIR}/include
    ${blasfeo_SOURCE_DIR}/include
    ${hpipm_SOURCE_DIR}/include
)
target_link_libraries(ctrlpp_qp_solvers PUBLIC
    ctrlpp::bench_problems
    argmin::argmin
    osqp::osqpstatic
    qpOASES
    daqpstat
    hpipm
    blasfeo
    proxsuite::proxsuite
    ctrlpp_qp_solvers_tinympc
)

# The allocation sensor replaces the C allocator's entry points with glibc-forwarding ones, so it is
# compiled into each driver directly and only where glibc provides them.
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(ctrlpp_qp_solvers_sensor alloc_sensor_glibc.cpp)
else()
    set(ctrlpp_qp_solvers_sensor alloc_sensor_absent.cpp)
endif()

foreach(cohort dense sparse)
    add_executable(qp_solvers_${cohort} ${cohort}_main.cpp ${ctrlpp_qp_solvers_sensor})
    target_link_libraries(qp_solvers_${cohort} PRIVATE ctrlpp_qp_solvers)
endforeach()
