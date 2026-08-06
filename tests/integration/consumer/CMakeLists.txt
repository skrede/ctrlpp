cmake_minimum_required(VERSION 3.25)
project(ctrlpp_integration_test)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CTRLPP_CONSUMER_BACKEND "" CACHE STRING "Backend component to request from the installed ctrlpp package")

if (CTRLPP_CONSUMER_BACKEND)
    find_package(ctrlpp CONFIG REQUIRED COMPONENTS ${CTRLPP_CONSUMER_BACKEND})

    add_executable(integration_test backend_${CTRLPP_CONSUMER_BACKEND}.cpp)
    target_link_libraries(integration_test PRIVATE ctrlpp::${CTRLPP_CONSUMER_BACKEND})
else ()
    find_package(ctrlpp CONFIG REQUIRED)

    add_executable(integration_test main.cpp)
    target_link_libraries(integration_test PRIVATE ctrlpp::ctrlpp)
endif ()

enable_testing()
add_test(NAME integration COMMAND integration_test)
