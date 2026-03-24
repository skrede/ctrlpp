#ifndef HPP_GUARD_CTRLPP_MODEL_DISCRETISE_H
#define HPP_GUARD_CTRLPP_MODEL_DISCRETISE_H

/// @brief Continuous-to-discrete state-space conversion (ZOH, Tustin, Euler, RK4).
///
/// @cite franklin2015 -- Franklin et al., "Feedback Control of Dynamic Systems", 2015
/// @cite astrom2006 -- Astrom & Hagglund, "Advanced PID Control", 2006

#include "ctrlpp/types.h"
#include "ctrlpp/model/state_space.h"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include <cstddef>

namespace ctrlpp
{

struct zoh
{
};

struct tustin
{
};

struct forward_euler
{
};

struct backward_euler
{
};

// zoh discretisation using Van Loan augmented matrix exponential method.
// Forms the augmented matrix M = [[A*dt, B*dt], [0, 0]], computes exp(M),
// then extracts Ad and Bd from the result. Cd = C, Dd = D.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY> discretise(zoh, const continuous_state_space<Scalar, NX, NU, NY>& sys, Scalar dt)
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);
    constexpr int aug = nx + nu;

    // Form augmented matrix
    Eigen::Matrix<Scalar, aug, aug> M;
    M.setZero();
    M.template block<nx, nx>(0, 0) = sys.A * dt;
    M.template block<nx, nu>(0, nx) = sys.B * dt;

    // Compute matrix exponential
    Eigen::Matrix<Scalar, aug, aug> expM = M.exp();

    // Extract discretised matrices
    Matrix<Scalar, NX, NX> Ad = expM.template block<nx, nx>(0, 0);
    Matrix<Scalar, NX, NU> Bd = expM.template block<nx, nu>(0, nx);

    return {Ad, Bd, sys.C, sys.D};
}

// Convenience wrapper matching the generic discretise<Method>(...) signature.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY> discretise(const continuous_state_space<Scalar, NX, NU, NY>& sys, Scalar dt, zoh = {})
{
    return discretise(zoh{}, sys, dt);
}

} // namespace ctrlpp

#endif
