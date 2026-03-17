if (CTRLPP_CMAKE_FETCH_DEPS)
    FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
        GIT_SHALLOW TRUE
        EXCLUDE_FROM_ALL
    )
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(Eigen3)
else ()
    find_package(Eigen3 3.4 CONFIG REQUIRED)
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

install(TARGETS ctrlpp_eigen EXPORT ctrlppTargets)

install(DIRECTORY ${INCLUDE_DIRECTORY}/ctrlpp/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ctrlpp
)
