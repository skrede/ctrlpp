# `sentinel` switches the allocation sensor's two halves, the Eigen sentinel define and the
# malloc-wrapping link options, and nothing else.

function(ctrlpp_nucleo_image target sentinel)
    add_executable(${target}
        nucleo_main.cpp
        usart3_console.cpp
        syscalls.cpp
        boot_handshake.cpp
        init_stubs.cpp
        build_id.cpp
        alloc_sensor.cpp
        operator_new.cpp
        family_report.cpp
        sbrk_ceiling.cpp
        cycle_counter.cpp
        posture.cpp
        reset_cause.cpp
        telemetry.cpp
        timing_run.cpp
        stack_watermark_board.cpp
        "${CMSIS_SYSTEM}"
        "${CMSIS_STARTUP}")

    target_include_directories(${target} SYSTEM PRIVATE
        "${CTRLPP_INCLUDE_DIR}"
        "${CTRLPP_SHARED_DIR}"
        "${eigen3_SOURCE_DIR}"
        "${CMSIS_CORE_INC}"
        "${CMSIS_DEV_INC}"
        $<$<BOOL:${CTRLPP_MCU_PREDICTIVE}>:${CTRLPP_MCU_ARGMIN_INC}>)

    # EIGEN_ALLOCA must be explicit on newlib: under strict -std=c++20 Eigen's alloca detection finds
    # nothing, so every kernel temporary would heap-allocate. FPv5 is scalar, so
    # EIGEN_MAX_ALIGN_BYTES=8 (natural double) removes the aligned-allocation path.
    target_compile_definitions(${target} PRIVATE
        STM32H753xx
        EIGEN_DONT_VECTORIZE
        EIGEN_MAX_ALIGN_BYTES=8
        EIGEN_ALLOCA=__builtin_alloca
        EIGEN_STACK_ALLOCATION_LIMIT=8192
        CTRLPP_NO_EXCEPTIONS
        CTRLPP_MCU_PREDICTIVE=$<BOOL:${CTRLPP_MCU_PREDICTIVE}>
        CTRLPP_MCU_BOOT_HANDSHAKE_MS=${CTRLPP_MCU_BOOT_HANDSHAKE_MS}
        CTRLPP_MCU_ALLOC_SENTINEL=$<BOOL:${sentinel}>)

    target_compile_options(${target} PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:${CTRLPP_MCU_OPTIMIZE_FLAG} -fno-exceptions -fno-rtti -Wno-psabi>
        -ffunction-sections -fdata-sections
        -fstack-usage)

    # No system-call stub library is linked, so a system call syscalls.cpp does not define is a link
    # error rather than a call that fails at run time.
    target_link_options(${target} PRIVATE
        ${CTRLPP_M7_LINK_ARCH}
        -T "${CMAKE_CURRENT_SOURCE_DIR}/stm32h753zi_flash.ld"
        -nostartfiles
        --specs=nano.specs
        # newlib-nano trims printf; the float-capable variant serves the verdict lines' %g and %e.
        -u _printf_float
        # Hashed over the image's own contents, so the fingerprint a host capture checks changes
        # whenever the image does.
        -Wl,--build-id=sha1
        -Wl,--gc-sections
        -Wl,-Map=${target}.map)

    # Routes the image's own calls to the C allocation family through the
    # counting definitions in alloc_sensor.cpp.
    if(sentinel)
        target_link_options(${target} PRIVATE
            -Wl,--wrap=malloc,--wrap=free,--wrap=calloc,--wrap=realloc)
    endif()

    # The tuned clock tree is raised inside a wrapped SystemInit rather than from
    # main, so an image built without the option is the reference image bit for
    # bit.
    if(CTRLPP_MCU_TUNED_POSTURE)
        target_sources(${target} PRIVATE clock_tree.cpp clock_tree_report.cpp)
        target_compile_definitions(${target} PRIVATE CTRLPP_MCU_TUNED_POSTURE)
        target_link_options(${target} PRIVATE -Wl,--wrap=SystemInit)
    endif()

    # Without the dependency an edited reserve in the linker script leaves the
    # old image in place, and the capture would describe a layout that was never
    # linked.
    set_target_properties(${target} PROPERTIES
        SUFFIX ".elf"
        LINK_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/stm32h753zi_flash.ld")

    # Report size at build time (must fit 2 MB flash / 128 KB DTCM).
    find_program(CTRLPP_ARM_SIZE "${CTRLPP_ARM_TOOLCHAIN_PREFIX}size")
    if(CTRLPP_ARM_SIZE)
        add_custom_command(TARGET ${target} POST_BUILD
            COMMAND "${CTRLPP_ARM_SIZE}" "$<TARGET_FILE:${target}>"
            VERBATIM)
    endif()
endfunction()
