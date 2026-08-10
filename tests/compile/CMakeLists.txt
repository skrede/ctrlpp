# Compile-time tests -- targets added as headers are implemented

add_executable(concept_check_test concept_check_test.cpp)
target_link_libraries(concept_check_test PRIVATE ctrlpp::ctrlpp)
target_compile_options(concept_check_test PRIVATE ${CTRLPP_WARNING_FLAGS})
add_test(NAME concept_check_test COMMAND concept_check_test)

add_executable(embedded_core_float embedded_core_float.cpp)
target_link_libraries(embedded_core_float PRIVATE ctrlpp::ctrlpp)
target_compile_options(embedded_core_float PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_FLOAT_TIER_STRICT_FLAGS})
add_test(NAME embedded_core_float COMMAND embedded_core_float)

add_executable(moving_horizon_factory_compile moving_horizon_factory_compile.cpp)
target_link_libraries(moving_horizon_factory_compile PRIVATE ctrlpp::ctrlpp)
target_compile_options(moving_horizon_factory_compile PRIVATE ${CTRLPP_WARNING_FLAGS})
add_test(NAME moving_horizon_factory_compile COMMAND moving_horizon_factory_compile)
