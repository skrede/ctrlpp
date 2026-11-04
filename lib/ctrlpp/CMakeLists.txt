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

# /bigobj travels with the target because the consumer's own instantiations reach MSVC's
# 65,279-section object ceiling: the lightest unit measured uses a third of it, and one with more
# than one scalar type exceeds it. ELF and Mach-O have no ceiling reachable at these counts.
target_compile_options(ctrlpp
    INTERFACE
    $<$<CXX_COMPILER_ID:MSVC>:/bigobj>
)

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

    # The export name defaults to the real target name, so without it an installed ctrlpp::osqp
    # does not exist and the component reads as absent; ALIAS targets are not exported.
    set_target_properties(ctrlpp_osqp PROPERTIES EXPORT_NAME osqp)
endif ()

if (CTRLPP_BUILD_NLOPT)
    add_library(ctrlpp_nlopt INTERFACE)
    add_library(ctrlpp::nlopt ALIAS ctrlpp_nlopt)

    target_link_libraries(ctrlpp_nlopt
        INTERFACE
        ctrlpp::ctrlpp
        NLopt::nlopt
    )

    target_compile_definitions(ctrlpp_nlopt
        INTERFACE
        CTRLPP_HAS_NLOPT=1
    )

    set_target_properties(ctrlpp_nlopt PROPERTIES EXPORT_NAME nlopt)
endif ()

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

    set_target_properties(ctrlpp_argmin PROPERTIES EXPORT_NAME argmin)
endif ()

# A fetched dependency that installs nothing belongs to no export set, so exporting ctrlpp around
# it is illegal; CMake would report that only at generate time, naming the sub-build target
# rather than the option a user can change.
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

    # The gate flags raised here are ordinary variables and not PARENT_SCOPE writes, because
    # configure_package_config_file substitutes them in this same directory scope.
    foreach (_ctrlpp_backend IN ITEMS OSQP NLOPT ARGMIN)
        if (CTRLPP_BUILD_${_ctrlpp_backend})
            string(TOLOWER ${_ctrlpp_backend} _ctrlpp_backend_target)
            install(TARGETS ctrlpp_${_ctrlpp_backend_target}
                EXPORT ctrlppTargets
            )
            set(CTRLPP_${_ctrlpp_backend}_INSTALLED TRUE)
        endif ()
    endforeach ()

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

    # SameMinorVersion because a pre-1.0 minor release may break the API, which major-version
    # compatibility would hide. From 1.0 the fit is CMake's semantic-version mode, not major-version
    # compatibility; it needs a CMake 4.4 producer, above the declared floor.
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
