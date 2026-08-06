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
    add_library(ctrlpp_osqp INTERFACE)
    add_library(ctrlpp::osqp ALIAS ctrlpp_osqp)

    target_link_libraries(ctrlpp_osqp
        INTERFACE
        ctrlpp::ctrlpp
        osqp::osqpstatic
    )

    target_compile_definitions(ctrlpp_osqp
        INTERFACE
        CTRLPP_HAS_OSQP=1
    )

    # The export namespace is prefixed to the target's export name, which defaults to the
    # real target name, so without this a consumer's link against ctrlpp::osqp finds nothing
    # and the component reads as absent even though it was built and installed. ALIAS targets
    # are not exported and cannot stand in.
    set_target_properties(ctrlpp_osqp PROPERTIES EXPORT_NAME osqp)
endif ()

# --- NLopt solver backend (optional) ---
if (CTRLPP_BUILD_NLOPT)
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

# Raised here rather than left to CMake. A fetched dependency that installs nothing is a real,
# non-imported target belonging to no export set, so the export around it is illegal; CMake
# detects that only at generate time and its message names the sub-build target instead of the
# option a user can change.
if (CTRLPP_INSTALL)
    set(_ctrlpp_declining "")
    get_property(_ctrlpp_facts GLOBAL PROPERTY CTRLPP_DEPENDENCY_INSTALL_FACTS)
    foreach (_ctrlpp_fact IN LISTS _ctrlpp_facts)
        string(REPLACE "|" ";" _ctrlpp_field "${_ctrlpp_fact}")
        list(GET _ctrlpp_field 0 _ctrlpp_name)
        list(GET _ctrlpp_field 1 _ctrlpp_option)
        list(GET _ctrlpp_field 2 _ctrlpp_fetched)
        list(GET _ctrlpp_field 3 _ctrlpp_installs)
        if (_ctrlpp_fetched AND NOT _ctrlpp_installs)
            if (_ctrlpp_option)
                list(APPEND _ctrlpp_declining "${_ctrlpp_name} (turn ${_ctrlpp_option} on)")
            else ()
                list(APPEND _ctrlpp_declining "${_ctrlpp_name}")
            endif ()
        endif ()
    endforeach ()

    if (_ctrlpp_declining)
        string(REPLACE ";" ", " _ctrlpp_declining "${_ctrlpp_declining}")
        message(FATAL_ERROR
            "ctrlpp cannot be installed. These dependencies were fetched into this build and "
            "install nothing, so they belong to no export set and ctrlpp's own export cannot "
            "reference them: ${_ctrlpp_declining}. Either supply a discoverable version of each "
            "one and configure with CTRLPP_CMAKE_FETCH_DEPS=OFF, or configure with "
            "CTRLPP_INSTALL=OFF if this build was never meant to ship."
        )
    endif ()

    install(TARGETS ctrlpp
        EXPORT ctrlppTargets
    )

    install(DIRECTORY ${INCLUDE_DIRECTORY}/ctrlpp/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ctrlpp
    )

    # The gate flags below are ordinary variables and not PARENT_SCOPE writes, because
    # configure_package_config_file substitutes them in this same directory scope.
    if (CTRLPP_BUILD_OSQP)
        install(TARGETS ctrlpp_osqp
            EXPORT ctrlppTargets
        )
        set(CTRLPP_OSQP_INSTALLED TRUE)
    endif ()

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

    # The compatibility argument is the whole of the version contract this package publishes.
    # Under this mode a 0.3.6 install satisfies a consumer who asked for 0.3.0 and refuses one
    # who asked for 0.1.0. Major-version compatibility would accept the 0.1.0 request, which
    # contradicts this project's own policy that pre-1.0 API breakage is free -- a consumer
    # would be handed a package built against a different API and told it matched.
    #
    # The accepted cost is on the far side of 1.0: this mode is stricter than semantic
    # versioning requires there, because a 1.1 install would refuse a request for 1.0 even
    # though nothing broke. At the 1.0 release the intended replacement is CMake's
    # semantic-version compatibility mode, which special-cases a zero major component and is
    # the exact semantics wanted -- NOT a return to major-version compatibility. It is not
    # written today because it needs a CMake 4.4 producer, above this project's declared floor.
    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/ctrlppConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMinorVersion
    )

    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/ctrlppConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/ctrlppConfigVersion.cmake"
        DESTINATION ${ctrlpp_INSTALL_CMAKEDIR}
    )
endif ()
