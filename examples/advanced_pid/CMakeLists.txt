add_executable(advanced_pid main.cpp)
target_link_libraries(advanced_pid PRIVATE ctrlpp::eigen)
target_compile_options(advanced_pid PRIVATE ${CPPCTRL_WARNING_FLAGS})
