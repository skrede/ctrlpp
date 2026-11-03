#ifndef HPP_GUARD_CTRLPP_CONTROL_PLACE_TYPES_H
#define HPP_GUARD_CTRLPP_CONTROL_PLACE_TYPES_H

namespace ctrlpp
{

// The first two are structural, so no change of data can lift them; the rest are conditions on the data, kept
// apart so a caller is never sent to change the one thing that cannot help.
enum class place_error
{
    multi_input_not_supported,
    multi_output_not_supported,
    non_finite_input,
    poles_not_conjugate_symmetric,
    // For place_observer this reports that (A, C) is unobservable.
    uncontrollable_pair,
};

}

#endif
