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

# /bigobj travels with the target rather than staying on ctrlpp's own test and
# example targets, because the cost it covers is incurred by whoever INSTANTIATES
# the templates, and that is the consumer. A header-only library that ships a
# ceiling its own tests need but its users do not get hands them a link failure
# whose cause is in a build file they cannot see.
#
# The margin is not generous. Counting COMDAT groups with g++ and
# per-function/per-data sections as a stand-in for the object format's own count,
# one intercept plus one continuous tuple plus one scalar type is roughly 22,500
# against the 65,279 ceiling -- 0.345 of it, and that is the LIGHTEST of the
# thirteen units measured. Three continuous instantiations in a single
# translation unit already sit where care_error_test sits, and more than one
# scalar type in that unit is past the ceiling outright. A consumer reaches this
# with ordinary use, not with abuse.
#
# It is MSVC-only because the ceiling is: the ELF and Mach-O section counts this
# would correspond to are not reachable at these instantiation counts, so there
# is nothing to widen elsewhere and no flag is imposed on those toolchains.
target_compile_options(ctrlpp
    INTERFACE
    $<$<CXX_COMPILER_ID:MSVC>:/bigobj>
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

# --- NLopt solver backend (optional) ---
if (CTRLPP_BUILD_NLOPT)
    FetchContent_Declare(
        nlopt
        GIT_REPOSITORY https://github.com/stevengj/nlopt.git
        GIT_TAG v2.10.1
        EXCLUDE_FROM_ALL
    )
    set(NLOPT_PYTHON OFF CACHE BOOL "" FORCE)
    set(NLOPT_OCTAVE OFF CACHE BOOL "" FORCE)
    set(NLOPT_GUILE OFF CACHE BOOL "" FORCE)
    set(NLOPT_TESTS OFF CACHE BOOL "" FORCE)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(nlopt)

    # Suppress warnings from NLopt headers (third-party code)
    get_target_property(_nlopt_inc nlopt INTERFACE_INCLUDE_DIRECTORIES)
    if (_nlopt_inc)
        set_target_properties(nlopt PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_nlopt_inc}")
    endif ()

    add_library(ctrlpp_nlopt INTERFACE)
    add_library(ctrlpp::nlopt ALIAS ctrlpp_nlopt)

    target_link_libraries(ctrlpp_nlopt
        INTERFACE
        ctrlpp::ctrlpp
        nlopt
    )

    target_compile_definitions(ctrlpp_nlopt
        INTERFACE
        CTRLPP_HAS_NLOPT=1
    )
endif ()

# --- Argmin solver backend (optional) ---
if (CTRLPP_BUILD_ARGMIN)
    if (CTRLPP_ARGMIN_SOURCE_DIR)
        FetchContent_Declare(
            argmin
            SOURCE_DIR "${CTRLPP_ARGMIN_SOURCE_DIR}"
            EXCLUDE_FROM_ALL
        )
    else ()
        FetchContent_Declare(
            argmin
            GIT_REPOSITORY https://github.com/skrede/argmin.git
            GIT_TAG ${CTRLPP_ARGMIN_GIT_TAG}
            EXCLUDE_FROM_ALL
        )
    endif ()
    set(ARGMIN_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(ARGMIN_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(ARGMIN_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(argmin)

    add_library(ctrlpp_argmin INTERFACE)
    add_library(ctrlpp::argmin ALIAS ctrlpp_argmin)

    target_link_libraries(ctrlpp_argmin
        INTERFACE
        ctrlpp::ctrlpp
        argmin::argmin
    )

    target_compile_definitions(ctrlpp_argmin
        INTERFACE
        CTRLPP_HAS_ARGMIN=1
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
