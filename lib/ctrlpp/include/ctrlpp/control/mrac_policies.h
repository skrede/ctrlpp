#ifndef HPP_GUARD_CTRLPP_CONTROL_MRAC_POLICIES_H
#define HPP_GUARD_CTRLPP_CONTROL_MRAC_POLICIES_H

/// @brief MRAC robustification policy tag types: dead-zone, sigma-modification, e-modification.
///
/// @cite slotine1991 -- Slotine & Li, "Applied Nonlinear Control", 1991, Ch. 8

namespace ctrlpp
{

struct no_robustification
{
    struct options_t {};
};

struct dead_zone
{
    template <typename Scalar>
    struct options_t
    {
        Scalar threshold{};
    };
};

struct sigma_modification
{
    template <typename Scalar>
    struct options_t
    {
        Scalar sigma{};
    };
};

struct e_modification
{
    template <typename Scalar>
    struct options_t
    {
        Scalar delta{};
    };
};

}

#endif
