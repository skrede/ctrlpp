#ifndef HPP_GUARD_CTRLPP_SYSID_SYSID_TYPES_H
#define HPP_GUARD_CTRLPP_SYSID_SYSID_TYPES_H

// Every enumerator except model_not_identified is an exact precondition of the routine's index arithmetic or
// data layout, not a quality judgment about the fit.
// @cite ljung1999 -- Ljung, "System Identification: Theory for the User", 1999, Ch. 4

namespace ctrlpp
{

enum class sysid_error
{
    record_length_mismatch,
    record_not_single_row,
    too_few_samples,
    non_finite_sample,
    block_rows_below_order,
    model_not_identified,
};

}

#endif
