add_subdirectory(ctrlpp)

if(CTRLPP_ENABLE_EIGEN_POLICY)
    add_subdirectory(ctrlpp-eigen)
endif()

if(CTRLPP_BUILD_IMPLOT)
    add_subdirectory(ctrlpp-implot)
endif()
