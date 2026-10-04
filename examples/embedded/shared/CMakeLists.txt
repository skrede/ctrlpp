add_executable(ctrlpp_generate_golden generate_golden.cpp)
target_link_libraries(ctrlpp_generate_golden PRIVATE ctrlpp::ctrlpp)
target_compile_options(ctrlpp_generate_golden PRIVATE ${CTRLPP_WARNING_FLAGS})

# The predictive family's reference needs the optional nonlinear-programming
# backend, so its generator exists only in a tree that has it.
if (TARGET ctrlpp::argmin)
    add_executable(ctrlpp_generate_predictive_golden generate_predictive_golden.cpp)
    target_link_libraries(ctrlpp_generate_predictive_golden PRIVATE ctrlpp::argmin)
    target_compile_options(ctrlpp_generate_predictive_golden PRIVATE ${CTRLPP_WARNING_FLAGS})
endif ()
