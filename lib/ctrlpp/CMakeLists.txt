set(INCLUDE_DIRECTORY include)

add_library(ctrlpp INTERFACE)
add_library(ctrlpp::ctrlpp ALIAS ctrlpp)

target_include_directories(ctrlpp
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${INCLUDE_DIRECTORY}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS ctrlpp
    EXPORT ctrlppTargets
)

install(DIRECTORY ${INCLUDE_DIRECTORY}/ctrlpp/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ctrlpp
)

# --- CMake package config ---
set(ctrlpp_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/ctrlpp)

install(EXPORT ctrlppTargets
    FILE ctrlppTargets.cmake
    NAMESPACE ctrlpp::
    DESTINATION ${ctrlpp_INSTALL_CMAKEDIR}
)

configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/cmake/ctrlppConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/ctrlppConfig.cmake"
    INSTALL_DESTINATION ${ctrlpp_INSTALL_CMAKEDIR}
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/ctrlppConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/ctrlppConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/ctrlppConfigVersion.cmake"
    DESTINATION ${ctrlpp_INSTALL_CMAKEDIR}
)
