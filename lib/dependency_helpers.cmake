# One recorded fact per acquired dependency, in a global property rather than a scoped
# variable because the decision that reads them back is taken in a directory that ran none
# of these sites. An empty install_option means the dependency has no setting that would
# place it in this build's prefix at all, so it can never belong to an export set.
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

# Eigen's are the only third-party headers this library puts in its own interface, and the
# warning set every target here compiles with draws roughly twelve thousand diagnostics out of
# them the moment they stop being treated as system headers -- a wall of noise no source change
# explains. An imported target carries that treatment by default and the interface target this
# build makes asks for it, so the outcome is asserted rather than left to a default a policy could
# flip. A target a parent project supplies is held to the same requirement.
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
