# Provisions the nonlinear-programming backend the predictive family needs, as
# the leg provisions CMSIS and Eigen: a pre-provisioned source root skips the
# fetch, and otherwise the backend is fetched for its headers only.
#
# The leg is its own project and cannot see the root project's cache, so the
# revision and the repository are read out of the files that declare them.
# Restating either here would make two pins that can drift apart.

set(CTRLPP_MCU_ARGMIN_DIR "" CACHE PATH
    "Pre-provisioned argmin source root; when set, skips the argmin FetchContent")

function(ctrlpp_mcu_read_declared out file pattern)
    file(STRINGS "${file}" _lines REGEX "${pattern}")
    list(GET _lines 0 _line)
    string(REGEX MATCH "${pattern}" _match "${_line}")
    if(NOT CMAKE_MATCH_1)
        message(FATAL_ERROR "ctrlpp_nucleo: no match for '${pattern}' in ${file}")
    endif()
    set(${out} "${CMAKE_MATCH_1}" PARENT_SCOPE)
endfunction()

if(CTRLPP_MCU_ARGMIN_DIR)
    set(argmin_SOURCE_DIR "${CTRLPP_MCU_ARGMIN_DIR}")
    message(STATUS "ctrlpp_nucleo: using provided argmin at ${argmin_SOURCE_DIR} (fetch skipped)")
else()
    set(_ctrlpp_root "${CMAKE_CURRENT_SOURCE_DIR}/../../..")
    ctrlpp_mcu_read_declared(_argmin_tag "${_ctrlpp_root}/CMakeLists.txt"
        "^set\\(CTRLPP_ARGMIN_GIT_TAG \"([^\"]+)\"")
    ctrlpp_mcu_read_declared(_argmin_repository "${_ctrlpp_root}/lib/dependencies.cmake"
        "GIT_REPOSITORY[ \t]+([^ \t]*/argmin\\.git)")
    message(STATUS "ctrlpp_nucleo: fetching argmin ${_argmin_tag} from ${_argmin_repository} (the root project's pin)")
    # A commit pin cannot be fetched shallow, so there is no GIT_SHALLOW here.
    FetchContent_Declare(argmin
        GIT_REPOSITORY "${_argmin_repository}"
        GIT_TAG "${_argmin_tag}"
        EXCLUDE_FROM_ALL
        SOURCE_SUBDIR this-directory-does-not-exist)
    FetchContent_MakeAvailable(argmin)
endif()

set(CTRLPP_MCU_ARGMIN_INC "${argmin_SOURCE_DIR}/lib/argmin/include")
