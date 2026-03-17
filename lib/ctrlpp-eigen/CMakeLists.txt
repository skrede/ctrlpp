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
    find_package(Eigen3 3.4 CONFIG REQUIRED)

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

add_library(ctrlpp_eigen INTERFACE)
add_library(ctrlpp::eigen ALIAS ctrlpp_eigen)

target_include_directories(ctrlpp_eigen
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${INCLUDE_DIRECTORY}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(ctrlpp_eigen
    INTERFACE
    ctrlpp::ctrlpp
    Eigen3::Eigen
)

if (NOT CTRLPP_CMAKE_FETCH_DEPS)
    install(TARGETS ctrlpp_eigen EXPORT ctrlppTargets)

    install(DIRECTORY ${INCLUDE_DIRECTORY}/ctrlpp/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ctrlpp
    )
endif ()
