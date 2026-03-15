add_executable(basic_pid main.cpp)
target_link_libraries(basic_pid PRIVATE ctrlpp::eigen)
target_compile_options(basic_pid PRIVATE ${CPPCTRL_WARNING_FLAGS})
