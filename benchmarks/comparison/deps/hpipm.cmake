include(FetchContent)

FetchContent_Declare(blasfeo
    GIT_REPOSITORY https://github.com/giaf/blasfeo.git
    GIT_TAG 9628393214623b29e5dd2a8784f7423edbc93e04
    EXCLUDE_FROM_ALL)

FetchContent_Declare(hpipm
    GIT_REPOSITORY https://github.com/giaf/hpipm.git
    GIT_TAG 16a8db59f595bda642f460aa0f62e308f4b76db3
    EXCLUDE_FROM_ALL)

# BLASFEO and HPIPM read one shared cache variable, TARGET, with different vocabularies, so it is
# forced to each project's own value immediately before that project is made available.
set(LA HIGH_PERFORMANCE CACHE STRING "" FORCE)
set(MF PANELMAJ CACHE STRING "" FORCE)
set(BLASFEO_TESTING OFF CACHE BOOL "" FORCE)
set(BLASFEO_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BLASFEO_BENCHMARKS OFF CACHE BOOL "" FORCE)
set(TARGET X64_INTEL_HASWELL CACHE STRING "" FORCE)
FetchContent_MakeAvailable(blasfeo)

set(HPIPM_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(TARGET AVX CACHE STRING "" FORCE)
FetchContent_MakeAvailable(hpipm)

# Read from each project's own compile flags, since the shared TARGET cache entry only holds the
# value forced last.
get_directory_property(blasfeo_c_flags DIRECTORY ${blasfeo_SOURCE_DIR} DEFINITION CMAKE_C_FLAGS)
get_directory_property(hpipm_c_flags DIRECTORY ${hpipm_SOURCE_DIR} DEFINITION CMAKE_C_FLAGS)
string(REGEX MATCH "-DTARGET_[A-Z0-9_]+" blasfeo_target "${blasfeo_c_flags}")
string(REGEX MATCH "-DTARGET_[A-Z0-9_]+" hpipm_target "${hpipm_c_flags}")
string(REPLACE "-DTARGET_" "" blasfeo_target "${blasfeo_target}")
string(REPLACE "-DTARGET_" "" hpipm_target "${hpipm_target}")
message(STATUS "BLASFEO effective TARGET=${blasfeo_target} LA=${LA} MF=${MF}")
message(STATUS "HPIPM effective TARGET=${hpipm_target}")
if(NOT blasfeo_target STREQUAL "X64_INTEL_HASWELL")
    message(FATAL_ERROR "BLASFEO compiles for target '${blasfeo_target}', expected X64_INTEL_HASWELL")
endif()
if(NOT hpipm_target STREQUAL "AVX")
    message(FATAL_ERROR "HPIPM compiles for target '${hpipm_target}', expected AVX (an unknown target falls back to GENERIC)")
endif()
