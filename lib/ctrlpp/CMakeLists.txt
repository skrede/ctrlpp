if (CTRLPP_CMAKE_FETCH_DEPS)
    FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
        GIT_SHALLOW TRUE
        EXCLUDE_FROM_ALL
        SOURCE_SUBDIR this-directory-does-not-exist
    )
    FetchContent_MakeAvailable(Eigen3)
    add_library(eigen_headers INTERFACE)
    add_library(Eigen3::Eigen ALIAS eigen_headers)
    target_include_directories(eigen_headers SYSTEM INTERFACE
        "${eigen3_SOURCE_DIR}"
    )
else ()
    find_package(Eigen3 CONFIG REQUIRED)
    if (Eigen3_VERSION VERSION_LESS "3.4")
        message(FATAL_ERROR "ctrlpp requires Eigen >= 3.4, found ${Eigen3_VERSION}")
    endif ()

    get_target_property(_eigen_aliased Eigen3::Eigen ALIASED_TARGET)
    if (_eigen_aliased)
        set(_eigen_target ${_eigen_aliased})
    else ()
        set(_eigen_target Eigen3::Eigen)
    endif ()
    get_target_property(_eigen_inc ${_eigen_target} INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(${_eigen_target} PROPERTIES
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_eigen_inc}"
    )
endif ()

set(INCLUDE_DIRECTORY include)

add_library(ctrlpp INTERFACE)
add_library(ctrlpp::ctrlpp ALIAS ctrlpp)

target_include_directories(ctrlpp
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${INCLUDE_DIRECTORY}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(ctrlpp
    INTERFACE
    Eigen3::Eigen
)

# --- OSQP solver backend (optional) ---
if (CTRLPP_BUILD_OSQP)
    FetchContent_Declare(
        osqp
        GIT_REPOSITORY https://github.com/osqp/osqp.git
        GIT_TAG v1.0.0
        EXCLUDE_FROM_ALL
    )
    set(OSQP_BUILD_DEMO OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(osqp)

    add_library(ctrlpp_osqp INTERFACE)
    add_library(ctrlpp::osqp ALIAS ctrlpp_osqp)

    target_link_libraries(ctrlpp_osqp
        INTERFACE
        ctrlpp::ctrlpp
        osqpstatic
    )

    target_compile_definitions(ctrlpp_osqp
        INTERFACE
        CTRLPP_HAS_OSQP=1
    )
endif ()

if (NOT CTRLPP_CMAKE_FETCH_DEPS)
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
endif ()
