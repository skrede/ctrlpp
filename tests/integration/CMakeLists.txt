include(ExternalProject)

set(_install_prefix "${CMAKE_BINARY_DIR}/integration-install")

ExternalProject_Add(ctrlpp_install
    SOURCE_DIR "${CMAKE_SOURCE_DIR}"
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${_install_prefix}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCPPCTRL_ENABLE_EIGEN_POLICY=${CPPCTRL_ENABLE_EIGEN_POLICY}
        -DCPPCTRL_CMAKE_FETCH_DEPS=${CPPCTRL_CMAKE_FETCH_DEPS}
        -DCPPCTRL_BUILD_TESTS=OFF
        -DCPPCTRL_BUILD_EXAMPLES=OFF
    BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR>
    INSTALL_COMMAND ${CMAKE_COMMAND} --install <BINARY_DIR>
    TEST_COMMAND ""
)

ExternalProject_Add(ctrlpp_consumer
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/consumer"
    CMAKE_ARGS
        -DCMAKE_PREFIX_PATH=${_install_prefix}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
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
