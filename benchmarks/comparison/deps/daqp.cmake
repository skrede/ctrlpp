include(FetchContent)

FetchContent_Declare(daqp
    GIT_REPOSITORY https://github.com/darnstrom/daqp.git
    GIT_TAG 29a8987d8b9a596c1ab6acf84c17d38b068b73bc
    EXCLUDE_FROM_ALL)

# FASTER_MATH defaults ON and adds value-changing floating-point flags no other solver gets; the
# architecture comes from the tree-wide flags rather than from NATIVE_OPTIMIZATION.
set(FASTER_MATH OFF CACHE BOOL "" FORCE)
set(NATIVE_OPTIMIZATION OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(daqp)

get_directory_property(daqp_c_flags DIRECTORY ${daqp_SOURCE_DIR} DEFINITION CMAKE_C_FLAGS)
if(daqp_c_flags MATCHES "associative-math|march=native")
    message(FATAL_ERROR "DAQP still compiles with FASTER_MATH or NATIVE_OPTIMIZATION flags: ${daqp_c_flags}")
endif()
message(STATUS "DAQP FASTER_MATH=${FASTER_MATH} NATIVE_OPTIMIZATION=${NATIVE_OPTIMIZATION}")
