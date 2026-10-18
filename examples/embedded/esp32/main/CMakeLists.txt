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

# The stack probe runs in the evidence image only, so only that image compiles
# it, one translation unit per dimension: the six discrete solves in one unit
# needed 11 GB to compile, and each image's copy was compiled in parallel with
# the other's. The pairs are the ones stack_probe_dimension.h declares. The
# framework first runs this file in script mode for its requirements alone,
# where nothing can be generated.
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

# Linked whole so the replaced operator new is the one in the image: as an
# archive member it would be pulled only if something had already left the
# operator undefined when the archive was scanned, and the C++ library's own
# definition would otherwise win.
idf_component_register(
    SRCS ${ctrlpp_esp32_sources} ${ctrlpp_esp32_probe_sources}
    INCLUDE_DIRS "."
    REQUIRES esp_driver_uart
    PRIV_REQUIRES esp_app_format
    WHOLE_ARCHIVE
)

# Everything one image's copy of this component needs. Both images share every
# definition and option below, and every source but the stack probe's, which
# the evidence image alone compiles; `sentinel` switches the allocation
# sensor -- the Eigen sentinel and the counting operators -- and with it the
# image kind, which the sources read to run the stack probe and the task
# watchdog in the evidence image only. The timing copy is also built with
# NDEBUG, so the code it times carries neither Eigen's checks nor the C
# assertions, as the other board's timing image does; the evidence copy keeps
# the framework's assertions and counts Eigen's checks through the sentinel.
#
# SYSTEM keeps the strict warning surface from tripping on third-party Eigen and
# the ctrlpp template depth. Eigen lives outside ESP-IDF's component graph and is
# not published as a managed component; the build expects a local checkout
# exposed via CTRLPP_EIGEN_DIR. CTRLPP_NO_EXCEPTIONS is belt-and-suspenders: the
# -fno-exceptions posture already drives CTRLPP_HAS_EXCEPTIONS to 0 via
# __cpp_exceptions, and defining it makes the no-throw intent explicit at the
# build boundary.
#
# XTENSA LITERAL POOLS, AND WHY THE STACK PROBE NEEDS -mtext-section-literals.
#
# Xtensa loads a constant with L32R, a PC-relative load reaching only backwards
# and only a bounded distance. By default the toolchain gathers a function's
# constants into a separate .literal section, and once a single function's code
# grows past that reach the linker refuses it outright:
#
#   dangerous relocation: l32r: literal target out of range
#
# The discrete acceptance path instantiated at state dimension eight in float is
# past it. That is not a warning about size, it is a hard link failure, and
# without this flag the probe cannot even be BUILT at the dimension it most needs
# to test. Placing literals inline in .text keeps every load within reach.
#
# The cost is slightly larger and slightly less cache-friendly code. It does not
# move stack frames, which is what this component is measuring.
# Guarded on the COMPILER, not on CMAKE_SYSTEM_PROCESSOR: ESP-IDF does not set
# that variable, so a guard on it silently never fires. The option is Xtensa-only
# and the RISC-V ESP32 parts would reject it.
function(ctrlpp_esp32_image lib sentinel)
    target_include_directories(${lib} SYSTEM PRIVATE
        "${CTRLPP_ROOT}/lib/ctrlpp/include"
        "${CTRLPP_ROOT}/examples/embedded/shared")
    if (DEFINED ENV{CTRLPP_EIGEN_DIR})
        target_include_directories(${lib} SYSTEM PRIVATE "$ENV{CTRLPP_EIGEN_DIR}")
    endif ()
    target_compile_definitions(${lib} PRIVATE CTRLPP_NO_EXCEPTIONS CTRLPP_MCU_ALLOC_SENTINEL=${sentinel})
    if (NOT sentinel)
        target_compile_definitions(${lib} PRIVATE NDEBUG)
    endif ()
    target_compile_features(${lib} PRIVATE cxx_std_20)
    if (CMAKE_C_COMPILER MATCHES "xtensa")
        target_compile_options(${lib} PRIVATE -mtext-section-literals)
    endif ()
endfunction()

ctrlpp_esp32_image(${COMPONENT_LIB} 1)

# The timing image's copy of this component. It is created in this directory so
# it inherits the framework's compile options exactly as the component does,
# and it links what the component links; the project links it into the timing
# executable in the component's place.
add_library(ctrlpp_esp32_timing_main STATIC ${ctrlpp_esp32_sources})
get_target_property(ctrlpp_esp32_main_links ${COMPONENT_LIB} LINK_LIBRARIES)
target_link_libraries(ctrlpp_esp32_timing_main PRIVATE ${ctrlpp_esp32_main_links})
target_include_directories(ctrlpp_esp32_timing_main PRIVATE "${CMAKE_CURRENT_LIST_DIR}")
ctrlpp_esp32_image(ctrlpp_esp32_timing_main 0)

# Nearly every source here parses Eigen and the ctrlpp surface, and the heaviest,
# the probe's dimension-eight solve, peaks at 3.4 GiB under the v6.0.2
# toolchain, so a build that schedules one per core runs out of memory on a
# machine with less than 3.4 GiB a core. The pool is sized from the machine's own
# memory and binds only where memory rather than the core count limits the
# build. Ninja alone honors it.
cmake_host_system_information(RESULT ctrlpp_esp32_memory_mib QUERY TOTAL_PHYSICAL_MEMORY)
math(EXPR ctrlpp_esp32_heavy_jobs "${ctrlpp_esp32_memory_mib} / 3482")
if (ctrlpp_esp32_heavy_jobs LESS 1)
    set(ctrlpp_esp32_heavy_jobs 1)
endif ()
set_property(GLOBAL APPEND PROPERTY JOB_POOLS ctrlpp_esp32_heavy=${ctrlpp_esp32_heavy_jobs})
set_property(TARGET ${COMPONENT_LIB} ctrlpp_esp32_timing_main PROPERTY JOB_POOL_COMPILE ctrlpp_esp32_heavy)
