cmake_minimum_required(VERSION 3.28)
project(ctrlpp_validation)

set(CMAKE_CXX_STANDARD 23)
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
add_subdirectory(cases/n4sid_ident)
