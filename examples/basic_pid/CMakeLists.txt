add_executable(ctrlpp_basic_pid main.cpp)
target_link_libraries(ctrlpp_basic_pid PRIVATE ctrlpp::eigen)
target_compile_options(ctrlpp_basic_pid PRIVATE ${CTRLPP_WARNING_FLAGS})
