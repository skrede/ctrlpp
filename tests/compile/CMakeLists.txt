# Compile-time tests -- targets added as headers are implemented

add_executable(concept_check_test concept_check_test.cpp)
target_link_libraries(concept_check_test PRIVATE ctrlpp::ctrlpp)
target_compile_options(concept_check_test PRIVATE ${CTRLPP_WARNING_FLAGS})
add_test(NAME concept_check_test COMMAND concept_check_test)
