cmake_minimum_required(VERSION 3.25)

# MSan fuzz builds require Clang with MSan-instrumented libc++.
# See LLVM docs for building libc++ with -fsanitize=memory.
option(CTRLPP_BUILD_FUZZ_MSAN "Build fuzz targets with MSan instead of ASan (separate build)" OFF)

function(make_fuzz_target TARGET_NAME)
    add_executable(${TARGET_NAME} "${CMAKE_CURRENT_LIST_DIR}/${TARGET_NAME}.cpp")
    target_link_libraries(${TARGET_NAME} PRIVATE ctrlpp::ctrlpp)
    target_compile_options(${TARGET_NAME} PRIVATE
        -fsanitize=fuzzer,address,undefined
        -fsanitize-recover=undefined
        -fno-omit-frame-pointer
    )
    target_link_options(${TARGET_NAME} PRIVATE
        -fsanitize=fuzzer,address,undefined
    )
endfunction()

function(make_fuzz_target_msan TARGET_NAME)
    set(MSAN_TARGET "${TARGET_NAME}_msan")
    add_executable(${MSAN_TARGET} "${CMAKE_CURRENT_LIST_DIR}/${TARGET_NAME}.cpp")
    target_link_libraries(${MSAN_TARGET} PRIVATE ctrlpp::ctrlpp)
    target_compile_options(${MSAN_TARGET} PRIVATE
        -fsanitize=fuzzer,memory
        -fsanitize-memory-track-origins=2
        -fno-omit-frame-pointer
    )
    target_link_options(${MSAN_TARGET} PRIVATE
        -fsanitize=fuzzer,memory
    )
endfunction()

make_fuzz_target(fuzz_dare)
make_fuzz_target(fuzz_kalman)
make_fuzz_target(fuzz_pid)

if (CTRLPP_BUILD_FUZZ_MSAN)
    make_fuzz_target_msan(fuzz_dare)
    make_fuzz_target_msan(fuzz_kalman)
    make_fuzz_target_msan(fuzz_pid)
endif ()

if (CTRLPP_BUILD_OSQP)
    add_executable(fuzz_mpc "${CMAKE_CURRENT_LIST_DIR}/fuzz_mpc.cpp")
    target_link_libraries(fuzz_mpc PRIVATE ctrlpp::osqp)
    target_compile_options(fuzz_mpc PRIVATE
        -fsanitize=fuzzer,address,undefined
        -fsanitize-recover=undefined
        -fno-omit-frame-pointer
    )
    target_link_options(fuzz_mpc PRIVATE
        -fsanitize=fuzzer,address,undefined
    )

    if (CTRLPP_BUILD_FUZZ_MSAN)
        add_executable(fuzz_mpc_msan "${CMAKE_CURRENT_LIST_DIR}/fuzz_mpc.cpp")
        target_link_libraries(fuzz_mpc_msan PRIVATE ctrlpp::osqp)
        target_compile_options(fuzz_mpc_msan PRIVATE
            -fsanitize=fuzzer,memory
            -fsanitize-memory-track-origins=2
            -fno-omit-frame-pointer
        )
        target_link_options(fuzz_mpc_msan PRIVATE
            -fsanitize=fuzzer,memory
        )
    endif ()
endif ()

if (CTRLPP_BUILD_NLOPT)
    add_executable(fuzz_nmpc "${CMAKE_CURRENT_LIST_DIR}/fuzz_nmpc.cpp")
    target_link_libraries(fuzz_nmpc PRIVATE ctrlpp::nlopt)
    target_compile_options(fuzz_nmpc PRIVATE
        -fsanitize=fuzzer,address,undefined
        -fsanitize-recover=undefined
        -fno-omit-frame-pointer
    )
    target_link_options(fuzz_nmpc PRIVATE
        -fsanitize=fuzzer,address,undefined
    )

    if (CTRLPP_BUILD_FUZZ_MSAN)
        add_executable(fuzz_nmpc_msan "${CMAKE_CURRENT_LIST_DIR}/fuzz_nmpc.cpp")
        target_link_libraries(fuzz_nmpc_msan PRIVATE ctrlpp::nlopt)
        target_compile_options(fuzz_nmpc_msan PRIVATE
            -fsanitize=fuzzer,memory
            -fsanitize-memory-track-origins=2
            -fno-omit-frame-pointer
        )
        target_link_options(fuzz_nmpc_msan PRIVATE
            -fsanitize=fuzzer,memory
        )
    endif ()
endif ()
