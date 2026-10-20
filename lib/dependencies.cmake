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

include(${CMAKE_CURRENT_LIST_DIR}/dependency_helpers.cmake)

# A dependency whose target a parent project already defines is the parent's: it is used as-is,
# with no find_package and no fetch, and no fact is recorded for it because whether it reaches an
# install prefix is the parent's decision rather than this build's.
#
# The version is a floor on the find path only, which changes one behavior: with fetching off,
# a system Eigen older than the floor no longer stops the configure but falls back to the pinned
# fetch. That fallback belongs to no export set, so an install asked for on top of it is refused
# by name by the readback in lib/ctrlpp/ rather than producing an unusable package.
#
# Deliberately not pinned to a commit like the two backends below: GIT_SHALLOW can fetch a tag
# but not an arbitrary commit on every host.
if (NOT TARGET Eigen3::Eigen)
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

# --- OSQP solver backend (optional) ---
if (CTRLPP_BUILD_OSQP AND NOT TARGET osqp::osqpstatic)
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

    # The fetched sub-build defines the bare name only -- its export() call writes a file and
    # creates no in-build alias -- while a found NLopt arrives namespaced. As with osqp above,
    # the namespaced name is aliased onto the bare one and never the reverse.
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

    # No install option is recorded: nlopt declares its install rules unconditionally and
    # EXCLUDE_FROM_ALL makes CMake ignore all of them, so nothing a user can set puts a fetched
    # nlopt in this build's prefix.
    _ctrlpp_record_dependency_install(nlopt "" "${nlopt_SOURCE_DIR}")
endif ()

# --- Argmin solver backend (optional) ---
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

    # No conditional alias: argmin declares argmin::argmin as a real alias in its own listfile
    # and its config package publishes the same spelling, so one link line serves both paths.
    #
    # No install option either, and unlike the two above that is not a consequence of
    # EXCLUDE_FROM_ALL: argmin gates its whole install block on being the top-level project,
    # which a fetched argmin never is. An installed argmin is therefore the only way this
    # backend can join an export set.
    _ctrlpp_record_dependency_install(argmin "" "${argmin_SOURCE_DIR}")
endif ()
