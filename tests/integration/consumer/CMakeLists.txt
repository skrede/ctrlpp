cmake_minimum_required(VERSION 3.25)
project(ctrlpp_integration_test)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ctrlpp CONFIG REQUIRED)

add_executable(integration_test_core main.cpp)
target_link_libraries(integration_test_core PRIVATE ctrlpp::ctrlpp)
target_compile_definitions(integration_test_core PRIVATE TEST_CORE_ONLY)

if(TARGET ctrlpp::eigen)
    add_executable(integration_test_eigen main.cpp)
    target_link_libraries(integration_test_eigen PRIVATE ctrlpp::eigen)
endif()

enable_testing()
add_test(NAME core COMMAND integration_test_core)
if(TARGET ctrlpp::eigen)
    add_test(NAME eigen COMMAND integration_test_eigen)
endif()
