add_executable(ctrlpp_advanced_pid main.cpp)
target_link_libraries(ctrlpp_advanced_pid PRIVATE ctrlpp::eigen)
target_compile_options(ctrlpp_advanced_pid PRIVATE ${CTRLPP_WARNING_FLAGS})
