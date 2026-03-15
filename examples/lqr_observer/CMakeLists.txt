add_executable(lqr_observer main.cpp)
target_link_libraries(lqr_observer PRIVATE ctrlpp::eigen)
target_compile_options(lqr_observer PRIVATE ${CPPCTRL_WARNING_FLAGS})
