# Standalone flashable NUCLEO-H753ZI image: CMSIS-only bare-metal (no HAL/RTOS),
# hand-written linker script, register-level USART3 console, and a bare superloop
# driving the shared double control kernel.
#
# This is its OWN project() below the examples/ include boundary: it consumes the
# ctrlpp headers via an include path only and is NEVER add_subdirectory'd into the
# host build. Configure it separately with the arm-none-eabi toolchain file:
#
#   cmake -S examples/embedded/nucleo_h753zi -B build-arm \
#     -DCMAKE_TOOLCHAIN_FILE=$PWD/examples/embedded/nucleo_h753zi/cmake/toolchain-arm-none-eabi.cmake \
#     -DCMAKE_BUILD_TYPE=Release
#   cmake --build build-arm

cmake_minimum_required(VERSION 3.28)
# C is needed for the CMSIS system_stm32h7xx.c (SystemInit); ASM for the gcc
# startup .s.
project(ctrlpp_nucleo_h753zi C CXX)

enable_language(ASM)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)

# ctrlpp headers, consumed header-only, plus the shared board-agnostic kernel.
set(CTRLPP_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../lib/ctrlpp/include")
set(CTRLPP_SHARED_DIR  "${CMAKE_CURRENT_SOURCE_DIR}/../shared")

# --- CMSIS-Core 6.x + STM32H7 device pack (headers + gcc startup + system) ---
# An integrator porting this onto their own board almost always already has a
# CMSIS core + device pack from their SDK; set CTRLPP_MCU_CMSIS_CORE_DIR and/or
# CTRLPP_MCU_CMSIS_DEVICE_DIR to those source roots to bypass the clone. Re-check
# the release tags before bumping them.
set(CTRLPP_MCU_CMSIS_CORE_DIR "" CACHE PATH
    "Pre-provisioned CMSIS-Core source root; when set, skips the cmsis_core FetchContent")
set(CTRLPP_MCU_CMSIS_DEVICE_DIR "" CACHE PATH
    "Pre-provisioned STM32H7 CMSIS device source root; when set, skips the cmsis_device_h7 FetchContent")

set(_cmsis_to_fetch "")
if(CTRLPP_MCU_CMSIS_CORE_DIR)
    set(cmsis_core_SOURCE_DIR "${CTRLPP_MCU_CMSIS_CORE_DIR}")
    message(STATUS "ctrlpp_nucleo: using provided CMSIS-Core at ${cmsis_core_SOURCE_DIR} (fetch skipped)")
else()
    FetchContent_Declare(cmsis_core
        GIT_REPOSITORY https://github.com/ARM-software/CMSIS_6.git
        GIT_TAG v6.3.0
        GIT_SHALLOW TRUE
        EXCLUDE_FROM_ALL
        SOURCE_SUBDIR this-directory-does-not-exist)
    list(APPEND _cmsis_to_fetch cmsis_core)
endif()
if(CTRLPP_MCU_CMSIS_DEVICE_DIR)
    set(cmsis_device_h7_SOURCE_DIR "${CTRLPP_MCU_CMSIS_DEVICE_DIR}")
    message(STATUS "ctrlpp_nucleo: using provided CMSIS device pack at ${cmsis_device_h7_SOURCE_DIR} (fetch skipped)")
else()
    FetchContent_Declare(cmsis_device_h7
        GIT_REPOSITORY https://github.com/STMicroelectronics/cmsis-device-h7.git
        GIT_TAG v1.10.7
        GIT_SHALLOW TRUE
        EXCLUDE_FROM_ALL
        SOURCE_SUBDIR this-directory-does-not-exist)
    list(APPEND _cmsis_to_fetch cmsis_device_h7)
endif()
if(_cmsis_to_fetch)
    FetchContent_MakeAvailable(${_cmsis_to_fetch})
endif()

set(CMSIS_CORE_INC "${cmsis_core_SOURCE_DIR}/CMSIS/Core/Include")
set(CMSIS_DEV_INC  "${cmsis_device_h7_SOURCE_DIR}/Include")
set(CMSIS_STARTUP  "${cmsis_device_h7_SOURCE_DIR}/Source/Templates/gcc/startup_stm32h753xx.s")
set(CMSIS_SYSTEM   "${cmsis_device_h7_SOURCE_DIR}/Source/Templates/system_stm32h7xx.c")

# Eigen, pinned to 3.4.0. Fetch is the default so a bare configure needs no setup;
# an integrator with their own Eigen sets CTRLPP_MCU_EIGEN_DIR to that source root.
set(CTRLPP_MCU_EIGEN_DIR "" CACHE PATH
    "Pre-provisioned Eigen source root; when set, skips the Eigen FetchContent")
if(CTRLPP_MCU_EIGEN_DIR)
    set(eigen3_SOURCE_DIR "${CTRLPP_MCU_EIGEN_DIR}")
    message(STATUS "ctrlpp_nucleo: using provided Eigen at ${eigen3_SOURCE_DIR} (fetch skipped)")
else()
    FetchContent_Declare(Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
        GIT_SHALLOW TRUE
        EXCLUDE_FROM_ALL
        SOURCE_SUBDIR this-directory-does-not-exist)
    FetchContent_MakeAvailable(Eigen3)
endif()

# Cortex-M7 hard-float arch flags. The toolchain file applies these at compile
# time; the gcc driver also needs them on the link line to select the correct
# newlib-nano multilib.
set(CTRLPP_M7_LINK_ARCH
    -mcpu=cortex-m7 -mfpu=fpv5-d16 -mfloat-abi=hard -mthumb)

add_executable(ctrlpp_nucleo_h753zi
    nucleo_main.cpp
    usart3_console.cpp
    init_stubs.cpp
    "${CMSIS_SYSTEM}"
    "${CMSIS_STARTUP}")

target_include_directories(ctrlpp_nucleo_h753zi SYSTEM PRIVATE
    "${CTRLPP_INCLUDE_DIR}"
    "${CTRLPP_SHARED_DIR}"
    "${eigen3_SOURCE_DIR}"
    "${CMSIS_CORE_INC}"
    "${CMSIS_DEV_INC}")

# EIGEN_ALLOCA must be explicit on newlib targets: at strict -std=c++20 (extensions
# OFF) Eigen's alloca auto-detection finds nothing, so every internal kernel
# temporary would heap-allocate instead of landing on the stack. FPv5 is scalar,
# so EIGEN_MAX_ALIGN_BYTES=8 (natural double) removes the aligned-allocation path.
target_compile_definitions(ctrlpp_nucleo_h753zi PRIVATE
    STM32H753xx
    EIGEN_DONT_VECTORIZE
    EIGEN_MAX_ALIGN_BYTES=8
    EIGEN_ALLOCA=__builtin_alloca
    EIGEN_STACK_ALLOCATION_LIMIT=8192
    CTRLPP_NO_EXCEPTIONS)

target_compile_options(ctrlpp_nucleo_h753zi PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-Os -fno-exceptions -fno-rtti -Wno-psabi>
    -ffunction-sections -fdata-sections
    -fstack-usage)

target_link_options(ctrlpp_nucleo_h753zi PRIVATE
    ${CTRLPP_M7_LINK_ARCH}
    -T "${CMAKE_CURRENT_SOURCE_DIR}/stm32h753zi_flash.ld"
    -nostartfiles
    --specs=nano.specs --specs=nosys.specs
    # newlib-nano trims printf: pull in the float-capable variant so the CSV
    # renders %f instead of a literal "f".
    -u _printf_float
    -Wl,--gc-sections
    -Wl,-Map=ctrlpp_nucleo_h753zi.map)

set_target_properties(ctrlpp_nucleo_h753zi PROPERTIES SUFFIX ".elf")

# Report size at build time (must fit 2 MB flash / 128 KB DTCM).
find_program(CTRLPP_ARM_SIZE "${CTRLPP_ARM_TOOLCHAIN_PREFIX}size")
if(CTRLPP_ARM_SIZE)
    add_custom_command(TARGET ctrlpp_nucleo_h753zi POST_BUILD
        COMMAND "${CTRLPP_ARM_SIZE}" "$<TARGET_FILE:ctrlpp_nucleo_h753zi>"
        VERBATIM)
endif()
