cmake_minimum_required(VERSION 3.28)
project(ctrlpp_validation)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 CONFIG REQUIRED)

# ctrlpp is header-only; create a lightweight interface target
add_library(ctrlpp_validation_lib INTERFACE)
add_library(ctrlpp::ctrlpp ALIAS ctrlpp_validation_lib)
target_include_directories(ctrlpp_validation_lib INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/../lib/ctrlpp/include
)
target_link_libraries(ctrlpp_validation_lib INTERFACE Eigen3::Eigen)

# This project is a separate CMake project that no standing build tree
# configures, and until now it compiled with no warning flags at all. That gap
# was not theoretical: it hid one translation unit that did not compile under
# the flags the rest of the tree uses, and two cases that discarded a typed
# fallible result inside their simulation loops -- a step that did nothing while
# the case believed it succeeded, in the exact harness whose whole purpose is to
# produce numbers a reference implementation is compared against.
#
# -Werror=unused-result is the one promoted to an error rather than left as a
# warning, matching the root project: these executables have no test framework
# and no assertions, so a discarded failure here produces a plausible-looking
# CSV rather than a diagnostic anyone reads.
#
# Attached to the interface target so every case inherits it, rather than
# repeated in nineteen leaf files where one could quietly drift.
target_compile_options(ctrlpp_validation_lib INTERFACE
    $<$<CXX_COMPILER_ID:MSVC>:/W4 /we4834>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:
        -Wall
        -Wextra
        -Werror=unused-result
    >
)

# Control module
add_subdirectory(cases/dare_solution)
add_subdirectory(cases/care_solution)
add_subdirectory(cases/lqr_step_response)
add_subdirectory(cases/lqr_finite_horizon)
add_subdirectory(cases/lqi_tracking)
add_subdirectory(cases/pid_linear_step)
add_subdirectory(cases/pole_placement)

# Model module
add_subdirectory(cases/c2d_zoh)
add_subdirectory(cases/analysis_poles)
add_subdirectory(cases/tf_ss_conversion)

# Estimation module
add_subdirectory(cases/kalman_tracking)
add_subdirectory(cases/luenberger_observer)

# DSP module
add_subdirectory(cases/butterworth_filter)
add_subdirectory(cases/fir_filter)

# Trajectory module
add_subdirectory(cases/cubic_spline_natural)

# Lie groups module
add_subdirectory(cases/so3_exp_log)

# System identification module
add_subdirectory(cases/batch_arx_ident)
add_subdirectory(cases/moesp_ident)
