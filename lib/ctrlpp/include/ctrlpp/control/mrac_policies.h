#ifndef HPP_GUARD_CTRLPP_CONTROL_MRAC_POLICIES_H
#define HPP_GUARD_CTRLPP_CONTROL_MRAC_POLICIES_H

// @cite slotine1991 -- Slotine & Li, "Applied Nonlinear Control", 1991, Ch. 8

namespace ctrlpp
{

struct no_robustification
{
    struct config
    {
    };
};

struct dead_zone
{
    template <typename Scalar>
    struct config
    {
        Scalar threshold{};
    };
};

struct sigma_modification
{
    template <typename Scalar>
    struct config
    {
        Scalar sigma{};
    };
};

struct e_modification
{
    template <typename Scalar>
    struct config
    {
        Scalar delta{};
    };
};

}

#endif
