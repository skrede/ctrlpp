include(FetchContent)

# Set here rather than at the root because the variable is directory-scoped and is read
# inside FetchContent_Declare: from here it reaches lib/ctrlpp/ by inheritance and does not
# reach tests/. That containment is the contract. tests/ declares RapidCheck outside the
# CTRLPP_CMAKE_FETCH_DEPS fork and then adds a subdirectory under ${rapidcheck_SOURCE_DIR},
# which only the fetch path defines, so a find-package-first policy reaching the test tree
# would leave that path with nothing to add.
if (CTRLPP_CMAKE_FETCH_DEPS)
    set(FETCHCONTENT_TRY_FIND_PACKAGE_MODE NEVER)
else ()
    set(FETCHCONTENT_TRY_FIND_PACKAGE_MODE ALWAYS)
endif ()

# One recorded fact per acquired dependency, in a global property rather than a scoped
# variable because the decision that reads them back is taken in a directory that ran none
# of these sites. An empty install_option means the dependency has no setting that would
# place it in this build's prefix at all, so it can never belong to an export set.
function(_ctrlpp_record_dependency_install name install_option source_dir)
    set(_fetched FALSE)
    set(_installs FALSE)
    if (source_dir)
        set(_fetched TRUE)
    endif ()
    if (install_option)
        if (${install_option})
            set(_installs TRUE)
        endif ()
    endif ()
    set_property(GLOBAL APPEND PROPERTY CTRLPP_DEPENDENCY_INSTALL_FACTS
        "${name}|${install_option}|${_fetched}|${_installs}"
    )
endfunction()

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
    _ctrlpp_record_dependency_install(Eigen3 "" "${eigen3_SOURCE_DIR}")
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

# --- OSQP solver backend (optional) ---
if (CTRLPP_BUILD_OSQP)
    # No version in FIND_PACKAGE_ARGS: osqp's own listfile declares no project version, so
    # every osqp config package reports 0.0.0 -- the distribution one and the one built from
    # the commit pinned here alike -- and any floor would make the find path fail for
    # everybody. GLOBAL is load-bearing: without it a found osqp arrives directory-scoped and
    # is invisible in the directories where the adapter's consumers define their targets.
    FetchContent_Declare(
        osqp
        GIT_REPOSITORY https://github.com/osqp/osqp.git
        GIT_TAG 236713ce9a56c182ac3230d52108f952afce1523  # v1.0.0
        EXCLUDE_FROM_ALL
        SYSTEM
        FIND_PACKAGE_ARGS NAMES osqp GLOBAL
    )
    set(OSQP_BUILD_DEMO OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(osqp)

    # A found osqp arrives namespaced from its own config package while the fetched sub-build
    # declares no alias at all, so one spelling has to serve both paths. The namespaced name
    # is aliased onto the bare one and never the reverse, because the namespaced spelling is
    # what gets written into the export file and what a consumer's find_dependency(osqp)
    # produces. Without this, the bare name is a legal plain library name in
    # target_link_libraries and reaches the linker as a raw flag rather than failing to
    # generate.
    if (NOT TARGET osqp::osqpstatic)
        add_library(osqp::osqpstatic ALIAS osqpstatic)
    endif ()

    # No install option is recorded: under EXCLUDE_FROM_ALL every install rule the fetched
    # sub-build declares is ignored, so nothing a user can set makes a fetched osqp reach
    # this build's prefix.
    _ctrlpp_record_dependency_install(osqp "" "${osqp_SOURCE_DIR}")
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
endif ()
