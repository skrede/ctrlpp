include(FetchContent)

# Built from source rather than found, so it compiles under the same flags as every other solver in
# the bench tree instead of the distribution's.
FetchContent_Declare(osqp
    GIT_REPOSITORY https://github.com/osqp/osqp.git
    GIT_TAG 236713ce9a56c182ac3230d52108f952afce1523
    EXCLUDE_FROM_ALL
    SYSTEM)

set(OSQP_BUILD_DEMO_EXE OFF CACHE BOOL "" FORCE)
set(OSQP_BUILD_UNITTESTS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(osqp)

if(NOT TARGET osqp::osqpstatic)
    add_library(osqp::osqpstatic ALIAS osqpstatic)
endif()
