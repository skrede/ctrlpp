# CMake toolchain file for bare-metal Cortex-M7 (NUCLEO-H753ZI) cross-compile.
#
# Used by the standalone nucleo_h753zi/ project() to build the full flashable
# image. It never touches the host build of the library or its tests. The
# compilers are not FetchContent-able; an arm-none-eabi GCC (14.x+) must be on
# PATH.

set(CMAKE_SYSTEM_NAME Generic)
set(CMAKE_SYSTEM_PROCESSOR arm)

# There is no startup/BSP at configure time, so a full hosted-executable link
# during CMake's compiler probe would fail; a static-library probe compiles
# without linking.
set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

set(CTRLPP_ARM_TOOLCHAIN_PREFIX "arm-none-eabi-" CACHE STRING
    "arm-none-eabi cross toolchain prefix")

find_program(CTRLPP_ARM_CC  "${CTRLPP_ARM_TOOLCHAIN_PREFIX}gcc" REQUIRED)
find_program(CTRLPP_ARM_CXX "${CTRLPP_ARM_TOOLCHAIN_PREFIX}g++" REQUIRED)

set(CMAKE_C_COMPILER   "${CTRLPP_ARM_CC}")
set(CMAKE_CXX_COMPILER "${CTRLPP_ARM_CXX}")
set(CMAKE_ASM_COMPILER "${CTRLPP_ARM_CC}")

# fpv5-d16 is the double+single-precision FPU the H753 actually has; the
# single-precision-only variant would be wrong for the double-precision leg.
set(CTRLPP_ARM_ARCH_FLAGS
    "-mcpu=cortex-m7 -mfpu=fpv5-d16 -mfloat-abi=hard -mthumb")

set(CMAKE_C_FLAGS_INIT   "${CTRLPP_ARM_ARCH_FLAGS}")
set(CMAKE_CXX_FLAGS_INIT "${CTRLPP_ARM_ARCH_FLAGS}")
set(CMAKE_ASM_FLAGS_INIT "${CTRLPP_ARM_ARCH_FLAGS}")

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)
