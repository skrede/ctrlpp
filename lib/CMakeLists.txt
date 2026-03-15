add_subdirectory(ctrlpp)

if(CPPCTRL_ENABLE_EIGEN_POLICY)
    add_subdirectory(ctrlpp-eigen)
endif()
