add_executable(ctrlpp_lqr_observer main.cpp)
target_link_libraries(ctrlpp_lqr_observer PRIVATE ctrlpp::eigen)
target_compile_options(ctrlpp_lqr_observer PRIVATE ${CTRLPP_WARNING_FLAGS})
