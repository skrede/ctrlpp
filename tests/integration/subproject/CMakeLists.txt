cmake_minimum_required(VERSION 3.28)
project(ctrlpp_subproject_test CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CTRLPP_SOURCE_DIR "" CACHE PATH "ctrlpp source tree to add as a subproject")
set(CTRLPP_FIXTURE_EIGEN_INCLUDE_DIR "" CACHE PATH "Eigen headers the parent's own Eigen3::Eigen carries")

# Hides every installed Eigen, so the parent's target below is the only Eigen ctrlpp can reach
# without fetching one, and a fetch is what the check after add_subdirectory catches.
set(CMAKE_DISABLE_FIND_PACKAGE_Eigen3 ON)

add_library(parent_eigen INTERFACE)
add_library(Eigen3::Eigen ALIAS parent_eigen)
target_include_directories(parent_eigen SYSTEM INTERFACE "${CTRLPP_FIXTURE_EIGEN_INCLUDE_DIR}")

add_subdirectory("${CTRLPP_SOURCE_DIR}" ctrlpp)

include(FetchContent)
FetchContent_GetProperties(Eigen3)
if (eigen3_POPULATED)
    message(FATAL_ERROR "ctrlpp fetched Eigen into ${eigen3_SOURCE_DIR} although its parent defines Eigen3::Eigen")
endif ()

add_executable(subproject_test "${CMAKE_CURRENT_SOURCE_DIR}/../consumer/main.cpp")
target_link_libraries(subproject_test PRIVATE ctrlpp::ctrlpp)

enable_testing()
add_test(NAME subproject COMMAND subproject_test)
