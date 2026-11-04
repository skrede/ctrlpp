#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_BOUNDS_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_BOUNDS_H

#include "ctrlpp/mpc/detail/nmpc_layout.h"

#include "ctrlpp/mpc/nmpc_config.h"

#include <Eigen/Dense>

#include <limits>
#include <cstddef>
#include <cstdint>

namespace ctrlpp::detail
{

template <typename Scalar>
struct nmpc_bounds
{
    Eigen::VectorX<Scalar> lower;
    Eigen::VectorX<Scalar> upper;
};

// x_0 is pinned by the initial-state equality alone: bounding it as well can make the NLP infeasible when x0
// sits on or over a bound, so the state bounds apply to stages 1..N.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC, std::size_t NTC>
inline void set_nmpc_state_bounds(const nmpc_config<Scalar, NX, NU, NC, NTC>& config, std::int32_t N, nmpc_bounds<Scalar>& bounds)
{
    constexpr int nx = static_cast<int>(NX);
    for(std::int32_t k = 1; k <= N; ++k)
    {
        if(config.x_min)
        {
            bounds.lower.segment(nmpc_layout::x_offset + k * nx, nx) = *config.x_min;
        }
        if(config.x_max)
        {
            bounds.upper.segment(nmpc_layout::x_offset + k * nx, nx) = *config.x_max;
        }
    }
}

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC, std::size_t NTC>
inline void set_nmpc_input_bounds(const nmpc_config<Scalar, NX, NU, NC, NTC>& config, std::int32_t N, Eigen::Index u_offset, nmpc_bounds<Scalar>& bounds)
{
    constexpr int nu = static_cast<int>(NU);
    for(std::int32_t k = 0; k < N; ++k)
    {
        if(config.u_min)
        {
            bounds.lower.segment(u_offset + k * nu, nu) = *config.u_min;
        }
        if(config.u_max)
        {
            bounds.upper.segment(u_offset + k * nu, nu) = *config.u_max;
        }
    }
}

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC, std::size_t NTC>
nmpc_bounds<Scalar> make_nmpc_variable_bounds(const nmpc_config<Scalar, NX, NU, NC, NTC>& config, const nmpc_layout& layout)
{
    nmpc_bounds<Scalar> bounds{Eigen::VectorX<Scalar>::Constant(layout.n_vars, -std::numeric_limits<Scalar>::infinity()),
            Eigen::VectorX<Scalar>::Constant(layout.n_vars, std::numeric_limits<Scalar>::infinity())};
    set_nmpc_state_bounds(config, layout.N, bounds);
    set_nmpc_input_bounds(config, layout.N, layout.u_offset, bounds);
    if(layout.has_path_slack)
    {
        bounds.lower.segment(layout.path_slack_offset, layout.n_path_slack).setZero();
    }
    if(layout.has_term_slack)
    {
        bounds.lower.segment(layout.term_slack_offset, layout.n_term_slack).setZero();
    }
    return bounds;
}

template <typename Scalar>
nmpc_bounds<Scalar> make_nmpc_constraint_bounds(const nmpc_layout& layout)
{
    nmpc_bounds<Scalar> bounds{Eigen::VectorX<Scalar>::Zero(layout.n_constraints), Eigen::VectorX<Scalar>::Zero(layout.n_constraints)};
    for(Eigen::Index i = layout.rate_start; i < layout.rate_start + layout.n_rate; ++i)
    {
        bounds.lower[i] = -std::numeric_limits<Scalar>::infinity();
    }
    for(Eigen::Index i = layout.path_con_start; i < layout.path_con_start + layout.n_path_con; ++i)
    {
        bounds.lower[i] = -std::numeric_limits<Scalar>::infinity();
    }
    for(Eigen::Index i = layout.term_con_start; i < layout.term_con_start + layout.n_term_con; ++i)
    {
        bounds.lower[i] = -std::numeric_limits<Scalar>::infinity();
    }
    return bounds;
}

}

#endif
