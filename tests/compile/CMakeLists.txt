# Compile-time tests -- targets added as headers are implemented

add_executable(concept_check_test concept_check_test.cpp)
target_link_libraries(concept_check_test PRIVATE ctrlpp::ctrlpp)
target_compile_options(concept_check_test PRIVATE ${CTRLPP_WARNING_FLAGS})
add_test(NAME concept_check_test COMMAND concept_check_test)

add_executable(embedded_core_float embedded_core_float.cpp)
target_link_libraries(embedded_core_float PRIVATE ctrlpp::ctrlpp)
target_compile_options(embedded_core_float PRIVATE ${CTRLPP_WARNING_FLAGS} ${CTRLPP_FLOAT_TIER_STRICT_FLAGS})
add_test(NAME embedded_core_float COMMAND embedded_core_float)

if (CTRLPP_BUILD_ARGMIN)
    add_executable(embedded_predictive_double embedded_predictive_double.cpp)
    target_link_libraries(embedded_predictive_double PRIVATE ctrlpp::argmin)
    target_compile_options(embedded_predictive_double PRIVATE ${CTRLPP_WARNING_FLAGS})
    add_test(NAME embedded_predictive_double COMMAND embedded_predictive_double)
endif ()

add_executable(moving_horizon_factory_compile moving_horizon_factory_compile.cpp)
target_link_libraries(moving_horizon_factory_compile PRIVATE ctrlpp::ctrlpp)
target_compile_options(moving_horizon_factory_compile PRIVATE ${CTRLPP_WARNING_FLAGS})
add_test(NAME moving_horizon_factory_compile COMMAND moving_horizon_factory_compile)

# The linear MPC's constructor is public so create can build it in place and set its solver up there; the passkey alone keeps an unvalidated, unset controller out.
ctrlpp_try_canary(CTRLPP_MPC_PASSKEY_BASELINE_COMPILED mpc_passkey_canary.cpp)
if (NOT CTRLPP_MPC_PASSKEY_BASELINE_COMPILED)
    message(FATAL_ERROR
        "The mpc passkey canary is not a valid positive control: "
        "its create arm does not compile.\n"
        "${CTRLPP_MPC_PASSKEY_BASELINE_COMPILED_OUTPUT}")
endif ()
ctrlpp_try_canary(CTRLPP_MPC_DIRECT_CONSTRUCTION_COMPILED mpc_passkey_canary.cpp -DCTRLPP_CANARY_MPC_DIRECT_CONSTRUCTION)
set(output "${CTRLPP_MPC_DIRECT_CONSTRUCTION_COMPILED_OUTPUT}")
if (CTRLPP_MPC_DIRECT_CONSTRUCTION_COMPILED)
    message(FATAL_ERROR "A direct mpc construction compiled; the passkey is not enforced.\n${output}")
elseif (NOT output MATCHES "validated_tag::validated_tag|constructor of class '[^']*::validated_tag'")
    message(FATAL_ERROR "A direct mpc construction failed for another reason than the passkey.\n${output}")
endif ()
