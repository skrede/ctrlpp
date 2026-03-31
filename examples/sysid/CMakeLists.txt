add_executable(rls_example rls_example.cpp)
target_link_libraries(rls_example PRIVATE ctrlpp::ctrlpp)
target_compile_options(rls_example PRIVATE ${CTRLPP_WARNING_FLAGS})

add_executable(batch_arx_example batch_arx_example.cpp)
target_link_libraries(batch_arx_example PRIVATE ctrlpp::ctrlpp)
target_compile_options(batch_arx_example PRIVATE ${CTRLPP_WARNING_FLAGS})
