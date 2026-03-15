#ifndef HPP_GUARD_CPPCTRL_EIGEN_DISCRETISE_H
#define HPP_GUARD_CPPCTRL_EIGEN_DISCRETISE_H

#include "ctrlpp/eigen_linalg.h"

#include <ctrlpp/discretise.h>
#include <ctrlpp/state_space.h>

#include <unsupported/Eigen/MatrixFunctions>

#include <cstddef>

namespace ctrlpp {

// ZOH discretisation using Van Loan augmented matrix exponential method.
// Forms the augmented matrix M = [[A*dt, B*dt], [0, 0]], computes exp(M),
// then extracts Ad and Bd from the result. Cd = C, Dd = D.
template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
auto discretise(ZOH, const ContinuousStateSpace<EigenLinalgPolicy, Scalar, NX, NU, NY>& sys,
                Scalar dt) -> DiscreteStateSpace<EigenLinalgPolicy, Scalar, NX, NU, NY>
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
    typename EigenLinalgPolicy::template matrix_type<Scalar, NX, NX> Ad =
        expM.template block<nx, nx>(0, 0);
    typename EigenLinalgPolicy::template matrix_type<Scalar, NX, NU> Bd =
        expM.template block<nx, nu>(0, nx);

    return {Ad, Bd, sys.C, sys.D};
}

// Convenience wrapper matching the generic discretise<Method>(...) signature.
template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
auto discretise(const ContinuousStateSpace<EigenLinalgPolicy, Scalar, NX, NU, NY>& sys,
                Scalar dt, ZOH = {}) -> DiscreteStateSpace<EigenLinalgPolicy, Scalar, NX, NU, NY>
{
    return discretise(ZOH{}, sys, dt);
}

}

#endif
