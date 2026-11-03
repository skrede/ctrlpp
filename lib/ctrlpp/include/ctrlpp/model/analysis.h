#ifndef HPP_GUARD_CTRLPP_MODEL_ANALYSIS_H
#define HPP_GUARD_CTRLPP_MODEL_ANALYSIS_H

// Each predicate answers whether the property is provable from the matrices it reads, so a NaN or an
// infinity in any of them is a negative answer; the finiteness guards are exact, with no tolerance.
// @cite anderson1990 -- Anderson & Moore, "Optimal Control: Linear Quadratic Methods", 1990
// @cite kailath1980 -- Kailath, "Linear Systems", 1980, Sec. 2.4 / Sec. 3.4 (Kalman rank tests)
// @cite franklin2015 -- Franklin et al., "Feedback Control of Dynamic Systems", 2015, Ch. 3 (stability criteria)

#include "ctrlpp/model/state_space.h"

#include "ctrlpp/model/detail/analysis_checks.h"

#include <Eigen/Eigenvalues>

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>

namespace ctrlpp
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
std::array<std::complex<Scalar>, NX> poles(const continuous_state_space<Scalar, NX, NU, NY>& sys)
{
    return detail::spectrum<Scalar, NX>(sys.A);
}

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
std::array<std::complex<Scalar>, NX> poles(const discrete_state_space<Scalar, NX, NU, NY>& sys)
{
    return detail::spectrum<Scalar, NX>(sys.A);
}

// @cite franklin2015 -- Franklin et al., "Feedback Control of Dynamic Systems", 2015, Ch. 3 (Hurwitz criterion)
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
bool is_stable(const continuous_state_space<Scalar, NX, NU, NY>& sys)
{
    if(!sys.A.allFinite())
        return false;

    auto p = poles(sys);
    for(const auto& pole : p)
    {
        if(!detail::is_determinate(pole))
            return false;
        if(pole.real() >= Scalar{0})
            return false;
    }
    return true;
}

// @cite astrom1997 -- Astrom & Wittenmark, "Computer-Controlled Systems", 3rd ed., 1997, Sec. 3.3 (unit-circle criterion)
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
bool is_stable(const discrete_state_space<Scalar, NX, NU, NY>& sys)
{
    if(!sys.A.allFinite())
        return false;

    auto p = poles(sys);
    for(const auto& pole : p)
    {
        if(!detail::is_determinate(pole))
            return false;
        if(std::abs(pole) >= Scalar{1})
            return false;
    }
    return true;
}

// @cite kailath1980 -- Kailath, "Linear Systems", 1980, Sec. 2.4 (Kalman rank test)
template <typename Scalar, std::size_t NX, std::size_t NU>
bool is_controllable(const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NX, NU>& B)
{
    if(!A.allFinite() || !B.allFinite())
        return false;

    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);

    Matrix<Scalar, NX, NX * NU> C;
    C.setZero();
    C.template block<nx, nu>(0, 0) = B;

    Matrix<Scalar, NX, NU> AkB = B;
    for(Eigen::Index k = 1; k < nx; ++k)
    {
        AkB = (A * AkB).eval();
        C.template block<nx, nu>(0, k * nu) = AkB;
    }

    return Eigen::FullPivLU<Matrix<Scalar, NX, NX * NU>>(C).rank() == Eigen::Index(NX);
}

// @cite kailath1980 -- Kailath, "Linear Systems", 1980, Sec. 3.4 (dual Kalman rank test)
template <typename Scalar, std::size_t NX, std::size_t NY>
bool is_observable(const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NY, NX>& C)
{
    if(!A.allFinite() || !C.allFinite())
        return false;

    constexpr int nx = static_cast<int>(NX);
    constexpr int ny = static_cast<int>(NY);

    Matrix<Scalar, NX * NY, NX> O;
    O.setZero();
    O.template block<ny, nx>(0, 0) = C;

    Matrix<Scalar, NY, NX> CAk = C;
    for(Eigen::Index k = 1; k < nx; ++k)
    {
        CAk = (CAk * A).eval();
        O.template block<ny, nx>(k * ny, 0) = CAk;
    }

    return Eigen::FullPivLU<Matrix<Scalar, NX * NY, NX>>(O).rank() == Eigen::Index(NX);
}

// The closed-loop matrix is checked as well as the operands: a finite triple can still overflow through
// B*K, and an infinity minus an infinity leaves a NaN behind.
template <typename Scalar, std::size_t NX, std::size_t NU>
bool is_stable_closed_loop(const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NX, NU>& B, const Matrix<Scalar, NU, NX>& K)
{
    if(!A.allFinite() || !B.allFinite() || !K.allFinite())
        return false;

    constexpr int n = static_cast<int>(NX);
    auto Acl = (A - B * K).eval();
    if(!Acl.allFinite())
        return false;

    Eigen::EigenSolver<Eigen::Matrix<Scalar, n, n>> solver(Acl, false);
    if(solver.info() != Eigen::Success)
        return false;

    auto evals = solver.eigenvalues();
    for(Eigen::Index i = 0; i < n; ++i)
        if(std::abs(evals(i)) >= Scalar{1})
            return false;
    return true;
}

// The observer matrix is checked as well as the operands, for the same reason as the closed-loop test.
template <typename Scalar, std::size_t NX, std::size_t NY>
bool is_stable_observer(const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NX, NY>& L, const Matrix<Scalar, NY, NX>& C)
{
    if(!A.allFinite() || !L.allFinite() || !C.allFinite())
        return false;

    constexpr int n = static_cast<int>(NX);
    auto Aobs = (A - L * C).eval();
    if(!Aobs.allFinite())
        return false;

    Eigen::EigenSolver<Eigen::Matrix<Scalar, n, n>> solver(Aobs, false);
    if(solver.info() != Eigen::Success)
        return false;

    auto evals = solver.eigenvalues();
    for(Eigen::Index i = 0; i < n; ++i)
        if(std::abs(evals(i)) >= Scalar{1})
            return false;
    return true;
}
}

#endif
