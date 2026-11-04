include(FetchContent)

# Set here, not at the root, so this directory-scoped policy reaches lib/ctrlpp/ but not tests/.
# tests/ adds RapidCheck from ${rapidcheck_SOURCE_DIR}, which only the fetch path defines.
if (CTRLPP_CMAKE_FETCH_DEPS)
    set(FETCHCONTENT_TRY_FIND_PACKAGE_MODE NEVER)
else ()
    set(FETCHCONTENT_TRY_FIND_PACKAGE_MODE ALWAYS)
endif ()

include(${CMAKE_CURRENT_LIST_DIR}/dependency_helpers.cmake)

# A parent project's Eigen target is used as-is and records no install fact: whether it reaches an
# install prefix is the parent's decision.
if (NOT TARGET Eigen3::Eigen)
    # The version is a floor on the find path only: an older system Eigen falls back to the fetch,
    # which an install then refuses by name. A tag rather than a commit pin, because GIT_SHALLOW
    # cannot fetch an arbitrary commit on every host.
    FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
        GIT_SHALLOW TRUE
        EXCLUDE_FROM_ALL
        SOURCE_SUBDIR this-directory-does-not-exist
        SYSTEM
        FIND_PACKAGE_ARGS 3.4 CONFIG NAMES Eigen3 GLOBAL
    )
    FetchContent_MakeAvailable(Eigen3)

    # SOURCE_SUBDIR names a directory that does not exist so the fetch adds no subdirectory and
    # builds none of Eigen's own tree; the target the headers travel on is made here instead.
    if (NOT TARGET Eigen3::Eigen)
        add_library(eigen_headers INTERFACE)
        add_library(Eigen3::Eigen ALIAS eigen_headers)
        target_include_directories(eigen_headers SYSTEM INTERFACE
            "${eigen3_SOURCE_DIR}"
        )
    endif ()
    _ctrlpp_record_dependency_install(Eigen3 "" "${eigen3_SOURCE_DIR}")
endif ()
_ctrlpp_require_system_eigen(Eigen3::Eigen)

if (CTRLPP_BUILD_OSQP AND NOT TARGET osqp::osqpstatic)
    # No version in FIND_PACKAGE_ARGS: osqp's listfile declares no project version, so every osqp
    # package reports 0.0.0. Without GLOBAL a found osqp is directory-scoped and invisible where
    # the adapter's consumers define their targets.
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

    # A found osqp is namespaced while the fetched sub-build declares no alias, so the namespaced
    # name, the one the export file and a consumer's find_dependency use, aliases the bare one.
    # Linking the bare name instead would, where no such target exists, pass it to the linker as a
    # raw flag rather than fail to generate.
    if (NOT TARGET osqp::osqpstatic)
        add_library(osqp::osqpstatic ALIAS osqpstatic)
    endif ()

    # No install option is recorded: EXCLUDE_FROM_ALL ignores every install rule of the fetched
    # sub-build, so nothing a user can set puts a fetched osqp in this build's prefix.
    _ctrlpp_record_dependency_install(osqp "" "${osqp_SOURCE_DIR}")
endif ()

if (CTRLPP_BUILD_NLOPT AND NOT TARGET NLopt::nlopt)
    # The find name is capitalized because the installed config file is NLoptConfig.cmake,
    # while the content name stays lowercase to match the sub-build's own project name.
    FetchContent_Declare(
        nlopt
        GIT_REPOSITORY https://github.com/stevengj/nlopt.git
        GIT_TAG 9e44e525370646def8152e73bb5c53a6531f6f7e  # v2.10.1
        EXCLUDE_FROM_ALL
        SYSTEM
        FIND_PACKAGE_ARGS 2.10 NAMES NLopt GLOBAL
    )
    set(NLOPT_PYTHON OFF CACHE BOOL "" FORCE)
    set(NLOPT_OCTAVE OFF CACHE BOOL "" FORCE)
    set(NLOPT_GUILE OFF CACHE BOOL "" FORCE)
    set(NLOPT_TESTS OFF CACHE BOOL "" FORCE)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(nlopt)

    # The fetched sub-build defines only the bare name (its export() creates no in-build alias)
    # while a found NLopt is namespaced; as with osqp, the namespaced name aliases the bare one.
    if (NOT TARGET NLopt::nlopt)
        add_library(NLopt::nlopt ALIAS nlopt)
    endif ()

    # Only a fetched, non-imported library reaches this: a found one is imported and its
    # include directories already carry system treatment.
    if (TARGET nlopt)
        get_target_property(_nlopt_imported nlopt IMPORTED)
        get_target_property(_nlopt_inc nlopt INTERFACE_INCLUDE_DIRECTORIES)
        if (NOT _nlopt_imported AND _nlopt_inc)
            set_target_properties(nlopt PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_nlopt_inc}")
        endif ()
    endif ()

    # No install option is recorded: EXCLUDE_FROM_ALL ignores nlopt's install rules, so nothing a
    # user can set puts a fetched nlopt in this build's prefix.
    _ctrlpp_record_dependency_install(nlopt "" "${nlopt_SOURCE_DIR}")
endif ()

if (CTRLPP_BUILD_ARGMIN AND NOT TARGET argmin::argmin)
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
            SYSTEM
            FIND_PACKAGE_ARGS ${CTRLPP_ARGMIN_MIN_VERSION} NAMES argmin GLOBAL
        )
    endif ()
    set(ARGMIN_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(ARGMIN_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(ARGMIN_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(argmin)

    # No alias: argmin's listfile and its config package both define argmin::argmin. No install
    # option either: argmin installs only as the top-level project, which a fetched one never is.
    _ctrlpp_record_dependency_install(argmin "" "${argmin_SOURCE_DIR}")
endif ()
