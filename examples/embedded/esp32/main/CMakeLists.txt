set(CTRLPP_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../../..")

set(ctrlpp_esp32_sources
    app_main.cpp
    telemetry.cpp
    posture.cpp
    stack_probe.cpp
    reset_reason.cpp
    family_report.cpp
    alloc_sensor.cpp
    operator_new.cpp
    cycle_counter.cpp
    timing_run.cpp)

# Linked whole so the replaced operator new is the one in the image: as an
# archive member it would be pulled only if something had already left the
# operator undefined when the archive was scanned, and the C++ library's own
# definition would otherwise win.
idf_component_register(
    SRCS ${ctrlpp_esp32_sources}
    INCLUDE_DIRS "."
    REQUIRES esp_driver_uart
    PRIV_REQUIRES esp_app_format
    WHOLE_ARCHIVE
)

# Everything one image's copy of this component needs. Both images share every
# source, definition and option below; `sentinel` switches the allocation
# sensor -- the Eigen sentinel and the counting operators -- and nothing else.
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
