# Compile-time tests -- targets added as headers are implemented

add_executable(concept_conformance_test concept_conformance_test.cpp)
target_link_libraries(concept_conformance_test PRIVATE ctrlpp::ctrlpp)
target_include_directories(concept_conformance_test PRIVATE ${CMAKE_SOURCE_DIR}/tests/common)
target_compile_options(concept_conformance_test PRIVATE ${CTRLPP_WARNING_FLAGS})

if (CTRLPP_ENABLE_EIGEN_POLICY)
    add_executable(eigen_conformance_test eigen_conformance_test.cpp)
    target_link_libraries(eigen_conformance_test PRIVATE ctrlpp::eigen)
    target_compile_options(eigen_conformance_test PRIVATE ${CTRLPP_WARNING_FLAGS})
endif ()
