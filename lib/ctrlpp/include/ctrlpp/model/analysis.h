#ifndef HPP_GUARD_CTRLPP_MODEL_ANALYSIS_H
#define HPP_GUARD_CTRLPP_MODEL_ANALYSIS_H

/// @brief State-space analysis: poles, controllability, observability.
///
/// @cite anderson1990 -- Anderson & Moore, "Optimal Control: Linear Quadratic Methods", 1990

#include "ctrlpp/model/state_space.h"

#include <Eigen/Eigenvalues>

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>

namespace ctrlpp
{

// poles: returns eigenvalues of the A matrix (system poles).
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
std::array<std::complex<Scalar>, NX> poles(const continuous_state_space<Scalar, NX, NU, NY>& sys)

{
    constexpr int n = static_cast<int>(NX);
    Eigen::EigenSolver<Eigen::Matrix<Scalar, n, n>> solver(sys.A, false);
    auto evals = solver.eigenvalues();
    std::array<std::complex<Scalar>, NX> result;
    for(int i = 0; i < n; ++i)
        result[static_cast<std::size_t>(i)] = evals(i);
    return result;
}

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
std::array<std::complex<Scalar>, NX> poles(const discrete_state_space<Scalar, NX, NU, NY>& sys)
{
    constexpr int n = static_cast<int>(NX);
    Eigen::EigenSolver<Eigen::Matrix<Scalar, n, n>> solver(sys.A, false);
    auto evals = solver.eigenvalues();

    std::array<std::complex<Scalar>, NX> result;
    for(int i = 0; i < n; ++i)
        result[static_cast<std::size_t>(i)] = evals(i);
    return result;
}

// is_stable: continuous system is stable iff all poles have negative real part.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
bool is_stable(const continuous_state_space<Scalar, NX, NU, NY>& sys)
{
    auto p = poles(sys);
    for(const auto& pole : p)
        if(pole.real() >= Scalar{0})
            return false;
    return true;
}

// is_stable: discrete system is stable iff all poles have magnitude < 1.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
bool is_stable(const discrete_state_space<Scalar, NX, NU, NY>& sys)
{
    auto p = poles(sys);
    for(const auto& pole : p)
        if(std::abs(pole) >= Scalar{1})
            return false;
    return true;
}

// is_controllable: checks rank of controllability matrix [B, AB, A^2 B, ..., A^{n-1} B].
// Returns true if rank equals NX (full state controllability).
template <typename Scalar, std::size_t NX, std::size_t NU>
bool is_controllable(const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NX, NU>& B)
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);

    Matrix<Scalar, NX, NX * NU> C;
    C.setZero();

    // Set first block: C[:, 0:NU] = B
    C.template block<nx, nu>(0, 0) = B;

    Matrix<Scalar, NX, NU> AkB = B;
    for(std::size_t k = 1; k < NX; ++k)
    {
        AkB = (A * AkB).eval();
        C.template block<nx, nu>(0, static_cast<int>(k * NU)) = AkB;
    }

    return static_cast<std::size_t>(Eigen::FullPivLU<Matrix<Scalar, NX, NX * NU>>(C).rank()) == NX;
}

// is_observable: checks rank of observability matrix [C; CA; CA^2; ...; CA^{n-1}].
// Returns true if rank equals NX (full state observability).
template <typename Scalar, std::size_t NX, std::size_t NY>
bool is_observable(const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NY, NX>& C)
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int ny = static_cast<int>(NY);

    Matrix<Scalar, NX * NY, NX> O;
    O.setZero();

    // Set first block: O[0:NY, :] = C
    O.template block<ny, nx>(0, 0) = C;

    Matrix<Scalar, NY, NX> CAk = C;
    for(std::size_t k = 1; k < NX; ++k)
    {
        CAk = (CAk * A).eval();
        O.template block<ny, nx>(static_cast<int>(k * NY), 0) = CAk;
    }

    return static_cast<std::size_t>(Eigen::FullPivLU<Matrix<Scalar, NX * NY, NX>>(O).rank()) == NX;
}

// is_stable_closed_loop: checks if all eigenvalues of (A - B*K) are inside the unit circle.
// For discrete-time closed-loop stability verification.
template <typename Scalar, std::size_t NX, std::size_t NU>
bool is_stable_closed_loop(const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NX, NU>& B, const Matrix<Scalar, NU, NX>& K)
{
    constexpr int n = static_cast<int>(NX);
    auto Acl = (A - B * K).eval();
    Eigen::EigenSolver<Eigen::Matrix<Scalar, n, n>> solver(Acl, false);
    auto evals = solver.eigenvalues();
    for(int i = 0; i < n; ++i)
        if(std::abs(evals(i)) >= Scalar{1})
            return false;
    return true;
}

// is_stable_observer: checks if all eigenvalues of (A - L*C) are inside the unit circle.
// For discrete-time observer stability verification.
template <typename Scalar, std::size_t NX, std::size_t NY>
bool is_stable_observer(const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NX, NY>& L, const Matrix<Scalar, NY, NX>& C)
{
    constexpr int n = static_cast<int>(NX);
    auto Aobs = (A - L * C).eval();
    Eigen::EigenSolver<Eigen::Matrix<Scalar, n, n>> solver(Aobs, false);
    auto evals = solver.eigenvalues();
    for(int i = 0; i < n; ++i)
        if(std::abs(evals(i)) >= Scalar{1})
            return false;
    return true;
}

} // namespace ctrlpp

#endif
