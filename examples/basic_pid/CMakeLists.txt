# Console variant -- always built when examples are enabled
add_executable(ctrlpp_basic_pid_console console/main.cpp)
target_link_libraries(ctrlpp_basic_pid_console PRIVATE ctrlpp::eigen)
target_compile_options(ctrlpp_basic_pid_console PRIVATE ${CTRLPP_WARNING_FLAGS})

# Interactive variant -- only when implot is available
if (CTRLPP_BUILD_IMPLOT)
    add_executable(ctrlpp_basic_pid interactive/main.cpp)
    target_link_libraries(ctrlpp_basic_pid PRIVATE ctrlpp::implot ctrlpp::eigen)
    target_compile_options(ctrlpp_basic_pid PRIVATE ${CTRLPP_WARNING_FLAGS})
endif ()
