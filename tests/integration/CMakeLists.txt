include(ExternalProject)

set(_toolchain_args "")
if (CMAKE_TOOLCHAIN_FILE)
    list(APPEND _toolchain_args "-DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}")
endif ()

# Configured at test time against the Eigen headers this tree already uses, so it needs neither an
# install nor a network.
add_test(
    NAME integration_subproject
    COMMAND ${CMAKE_CTEST_COMMAND}
        --build-and-test "${CMAKE_CURRENT_SOURCE_DIR}/subproject" "${CMAKE_CURRENT_BINARY_DIR}/subproject"
        --build-generator "${CMAKE_GENERATOR}"
        --build-config $<CONFIG>
        --build-options
            "-DCTRLPP_SOURCE_DIR=${CMAKE_SOURCE_DIR}"
            "-DCTRLPP_FIXTURE_EIGEN_INCLUDE_DIR=$<TARGET_PROPERTY:Eigen3::Eigen,INTERFACE_INCLUDE_DIRECTORIES>"
            "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
            "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
            ${_toolchain_args}
        --test-command ${CMAKE_CTEST_COMMAND} --output-on-failure
)
set_tests_properties(integration_subproject PROPERTIES
    LABELS "integration"
)

# The install round trip below needs every dependency discoverable, which a fetched one is not.
if (CTRLPP_CMAKE_FETCH_DEPS)
    return()
endif ()

set(_install_prefix "${CMAKE_BINARY_DIR}/integration-install")

ExternalProject_Add(ctrlpp_install
    SOURCE_DIR "${CMAKE_SOURCE_DIR}"
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${_install_prefix}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCTRLPP_CMAKE_FETCH_DEPS=OFF
        -DCTRLPP_BUILD_TESTS=OFF
        -DCTRLPP_BUILD_EXAMPLES=OFF
        ${_toolchain_args}
    BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR>
    INSTALL_COMMAND ${CMAKE_COMMAND} --install <BINARY_DIR>
    TEST_COMMAND ""
)

ExternalProject_Add(ctrlpp_consumer
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/consumer"
    CMAKE_ARGS
        -DCMAKE_PREFIX_PATH=${_install_prefix}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        ${_toolchain_args}
    DEPENDS ctrlpp_install
    INSTALL_COMMAND ""
    TEST_COMMAND ""
)

ExternalProject_Get_Property(ctrlpp_consumer BINARY_DIR)

add_test(
    NAME integration_find_package
    COMMAND ${CMAKE_COMMAND} --build ${BINARY_DIR} --target test
)
set_tests_properties(integration_find_package PROPERTIES
    LABELS "integration"
    DEPENDS ctrlpp_consumer
)
