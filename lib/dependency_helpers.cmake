# A global property rather than a scoped variable, because the install decision that reads these
# facts runs in a directory none of the recording sites ran in. An empty install_option means no
# setting places the dependency in this build's prefix, so it can never join an export set.
function(_ctrlpp_record_dependency_install name install_option source_dir)
    set(_fetched FALSE)
    set(_installs FALSE)
    if (source_dir)
        set(_fetched TRUE)
    endif ()
    if (install_option)
        if (${install_option})
            set(_installs TRUE)
        endif ()
    endif ()
    set_property(GLOBAL APPEND PROPERTY CTRLPP_DEPENDENCY_INSTALL_FACTS
        "${name}|${install_option}|${_fetched}|${_installs}"
    )
endfunction()

# Eigen's are the only third-party headers in this library's interface, and without system-header
# treatment the project's warning set draws about twelve thousand diagnostics from them. The
# treatment is asserted, for a parent project's target too, rather than left to a default a policy
# could flip.
function(_ctrlpp_require_system_eigen target)
    get_target_property(_aliased ${target} ALIASED_TARGET)
    if (_aliased)
        set(target ${_aliased})
    endif ()
    get_target_property(_imported ${target} IMPORTED)
    get_target_property(_system_inc ${target} INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
    set(_no_system FALSE)
    if (_imported)
        get_target_property(_no_system ${target} IMPORTED_NO_SYSTEM)
    endif ()
    if (_no_system OR NOT (_imported OR _system_inc))
        message(FATAL_ERROR
            "the Eigen target ${target} would hand its headers to every consumer as ordinary "
            "includes instead of system includes, which buries this project's own diagnostics under "
            "third-party ones."
        )
    endif ()
endfunction()
