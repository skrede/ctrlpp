set(CTRLPP_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../../..")

idf_component_register(
    SRCS "app_main.cpp"
    INCLUDE_DIRS "."
    REQUIRES esp_driver_uart
)

# ctrlpp's public include tree plus the shared board-leg kernel. SYSTEM keeps the
# strict warning surface from tripping on third-party Eigen and the ctrlpp
# template depth.
target_include_directories(${COMPONENT_LIB} SYSTEM PRIVATE
    "${CTRLPP_ROOT}/lib/ctrlpp/include"
    "${CTRLPP_ROOT}/examples/embedded/shared")

# Eigen lives outside ESP-IDF's component graph and is not published as a
# managed component; the build expects a local checkout exposed via
# CTRLPP_EIGEN_DIR.
if (DEFINED ENV{CTRLPP_EIGEN_DIR})
    target_include_directories(${COMPONENT_LIB} SYSTEM PRIVATE "$ENV{CTRLPP_EIGEN_DIR}")
endif ()

# Belt-and-suspenders: the -fno-exceptions posture already drives
# CTRLPP_HAS_EXCEPTIONS to 0 via __cpp_exceptions, and defining this makes the
# no-throw intent explicit at the build boundary.
target_compile_definitions(${COMPONENT_LIB} PRIVATE CTRLPP_NO_EXCEPTIONS)

# C++20 baseline matches the ctrlpp project floor.
target_compile_features(${COMPONENT_LIB} PRIVATE cxx_std_20)

# XTENSA LITERAL POOLS, AND WHY THE STACK PROBE NEEDS THIS.
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
if (CMAKE_C_COMPILER MATCHES "xtensa")
    target_compile_options(${COMPONENT_LIB} PRIVATE -mtext-section-literals)
endif ()
