include(FetchContent)

FetchContent_Declare(qpoases
    GIT_REPOSITORY https://github.com/coin-or/qpOASES.git
    GIT_TAG 680f18c8ef0018a120e1604b769f056e8368df97
    EXCLUDE_FROM_ALL)

set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(QPOASES_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

# qpOASES ships cmake_minimum_required(VERSION 2.6); lift to 3.5 so the
# configure step does not error on modern CMake.
set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "" FORCE)

FetchContent_MakeAvailable(qpoases)

# qpOASES replaces CMAKE_CXX_FLAGS_RELEASE with its debug flags plus -O3, which drops the
# -DNDEBUG every other solver in the tree compiles with.
target_compile_definitions(qpOASES PRIVATE NDEBUG)
