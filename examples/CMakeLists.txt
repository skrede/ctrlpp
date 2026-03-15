# ctrlpp examples
# Add example subdirectories here as they are created.

if (CPPCTRL_ENABLE_EIGEN_POLICY)
    add_subdirectory(basic_pid)
    add_subdirectory(advanced_pid)
    add_subdirectory(lqr_observer)
endif ()
