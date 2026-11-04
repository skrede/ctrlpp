set(CTRLPP_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../../..")

set(ctrlpp_esp32_sources
    app_main.cpp
    telemetry.cpp
    posture.cpp
    reset_reason.cpp
    family_report.cpp
    alloc_sensor.cpp
    alloc_canary.cpp
    operator_new.cpp
    cycle_counter.cpp
    timing_run.cpp)

# One translation unit per dimension pair of stack_probe_dimension.h, since the six discrete solves
# in one unit needed 11 GB to compile. The framework first runs this file in script mode for its
# requirements alone, where nothing can be generated.
set(ctrlpp_esp32_probe_sources stack_probe.cpp)
if (NOT CMAKE_BUILD_EARLY_EXPANSION)
    foreach (dimension IN ITEMS 2x1 3x1 4x2 5x1 6x2 8x2)
        string(REPLACE "x" ";" pair ${dimension})
        list(GET pair 0 PROBE_NX)
        list(GET pair 1 PROBE_NU)
        set(unit "${CMAKE_CURRENT_BINARY_DIR}/stack_probe_nx${PROBE_NX}_nu${PROBE_NU}.cpp")
        configure_file(stack_probe_dimension.cpp.in "${unit}" @ONLY)
        list(APPEND ctrlpp_esp32_probe_sources "${unit}")
    endforeach ()
endif ()

# Linked whole so the replaced operator new is the image's: as an archive member it would be pulled
# only if the operator were already undefined, and the C++ library's own definition would win.
idf_component_register(
    SRCS ${ctrlpp_esp32_sources} ${ctrlpp_esp32_probe_sources}
    INCLUDE_DIRS "."
    REQUIRES esp_driver_uart
    PRIV_REQUIRES esp_app_format
    WHOLE_ARCHIVE
)

# `sentinel` switches the allocation sensor (the Eigen sentinel and the counting operators) and with
# it the image kind the sources read. The timing copy also gets NDEBUG, so the code it times carries
# neither Eigen's checks nor the C assertions.
function(ctrlpp_esp32_image lib sentinel)
    # SYSTEM keeps the strict warning surface from tripping on Eigen and the ctrlpp template depth.
    target_include_directories(${lib} SYSTEM PRIVATE
        "${CTRLPP_ROOT}/lib/ctrlpp/include"
        "${CTRLPP_ROOT}/examples/embedded/shared")
    # Eigen is outside ESP-IDF's component graph, so the build expects a local checkout exposed via
    # CTRLPP_EIGEN_DIR.
    if (DEFINED ENV{CTRLPP_EIGEN_DIR})
        target_include_directories(${lib} SYSTEM PRIVATE "$ENV{CTRLPP_EIGEN_DIR}")
    endif ()
    # -fno-exceptions already drives CTRLPP_HAS_EXCEPTIONS to 0; CTRLPP_NO_EXCEPTIONS states the
    # intent at the build boundary.
    target_compile_definitions(${lib} PRIVATE CTRLPP_NO_EXCEPTIONS CTRLPP_MCU_ALLOC_SENTINEL=${sentinel})
    if (NOT sentinel)
        target_compile_definitions(${lib} PRIVATE NDEBUG)
    endif ()
    target_compile_features(${lib} PRIVATE cxx_std_20)
    # Xtensa's L32R reaches a literal only a bounded distance backwards, and the float discrete
    # acceptance path at state dimension eight outgrows it ("dangerous relocation: l32r: literal target
    # out of range"); literals inline in .text keep every load in reach and move no stack frame.
    # Guarded on the compiler because ESP-IDF does not set CMAKE_SYSTEM_PROCESSOR and RISC-V parts
    # reject the option.
    if (CMAKE_C_COMPILER MATCHES "xtensa")
        target_compile_options(${lib} PRIVATE -mtext-section-literals)
    endif ()
endfunction()

ctrlpp_esp32_image(${COMPONENT_LIB} 1)

# The timing image's copy of this component, created here so that it inherits the framework's
# compile options exactly as the component does.
add_library(ctrlpp_esp32_timing_main STATIC ${ctrlpp_esp32_sources})
get_target_property(ctrlpp_esp32_main_links ${COMPONENT_LIB} LINK_LIBRARIES)
target_link_libraries(ctrlpp_esp32_timing_main PRIVATE ${ctrlpp_esp32_main_links})
target_include_directories(ctrlpp_esp32_timing_main PRIVATE "${CMAKE_CURRENT_LIST_DIR}")
ctrlpp_esp32_image(ctrlpp_esp32_timing_main 0)

# The heaviest unit, the probe's dimension-eight solve, peaks at 3.4 GiB under the v6.0.2
# toolchain, so the pool is sized from the machine's memory and binds where memory rather than the
# core count limits the build. Ninja alone honors it.
cmake_host_system_information(RESULT ctrlpp_esp32_memory_mib QUERY TOTAL_PHYSICAL_MEMORY)
math(EXPR ctrlpp_esp32_heavy_jobs "${ctrlpp_esp32_memory_mib} / 3482")
if (ctrlpp_esp32_heavy_jobs LESS 1)
    set(ctrlpp_esp32_heavy_jobs 1)
endif ()
set_property(GLOBAL APPEND PROPERTY JOB_POOLS ctrlpp_esp32_heavy=${ctrlpp_esp32_heavy_jobs})
set_property(TARGET ${COMPONENT_LIB} ctrlpp_esp32_timing_main PROPERTY JOB_POOL_COMPILE ctrlpp_esp32_heavy)
