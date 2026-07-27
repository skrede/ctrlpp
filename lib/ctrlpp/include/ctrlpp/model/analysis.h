#ifndef HPP_GUARD_CTRLPP_MODEL_ANALYSIS_H
#define HPP_GUARD_CTRLPP_MODEL_ANALYSIS_H

/// @brief State-space analysis: poles, controllability, observability.
///
/// Every predicate on this page answers whether the property is PROVABLE from
/// the matrices it is given, so an indeterminate input is a negative answer: a
/// matrix carrying a NaN or an infinity establishes nothing, and the predicate
/// returns false rather than reporting a property it cannot demonstrate. The
/// guards are exact finiteness tests on the matrices each function reads, with
/// no tolerance and no threshold.
///
/// @cite anderson1990 -- Anderson & Moore, "Optimal Control: Linear Quadratic Methods", 1990
/// @cite kailath1980 -- Kailath, "Linear Systems", 1980, Sec. 2.4 / Sec. 3.4 (Kalman rank tests)
/// @cite franklin2015 -- Franklin et al., "Feedback Control of Dynamic Systems", 2015, Ch. 3 (stability criteria)

#include "ctrlpp/model/state_space.h"

#include <Eigen/Eigenvalues>

#include <array>
#include <cmath>
#include <limits>
#include <complex>
#include <cstddef>

namespace ctrlpp
{

namespace detail
{

/// @brief Eigenvalues of a square matrix, or an all-NaN array when no spectrum
/// is computable from it.
///
/// Two conditions leave the eigenvalue vector without meaning: a non-finite
/// entry in the matrix, for which no eigenvalue is defined, and a solver that
/// reports anything other than `Eigen::Success`, which leaves its output
/// unspecified. Both fill the whole result with a quiet NaN, so a caller cannot
/// mistake an unspecified or partially populated vector for a computed
/// spectrum. There is no tolerance here: the conditions are exact.
template <typename Scalar, std::size_t NX>
std::array<std::complex<Scalar>, NX> spectrum(const Matrix<Scalar, NX, NX>& A)
{
    constexpr int n = static_cast<int>(NX);
    constexpr Scalar nan = std::numeric_limits<Scalar>::quiet_NaN();

    std::array<std::complex<Scalar>, NX> result;
    result.fill(std::complex<Scalar>{nan, nan});

    if(!A.allFinite())
        return result;

    Eigen::EigenSolver<Eigen::Matrix<Scalar, n, n>> solver(A, false);
    if(solver.info() != Eigen::Success)
        return result;

    auto evals = solver.eigenvalues();
    for(int i = 0; i < n; ++i)
        result[static_cast<std::size_t>(i)] = evals(i);
    return result;
}

/// @brief True when a complex pole is a computed value rather than the NaN
/// marker `spectrum` emits for an indeterminate input.
template <typename Scalar>
bool is_determinate(const std::complex<Scalar>& pole)
{
    return std::isfinite(pole.real()) && std::isfinite(pole.imag());
}

}

/// @brief Returns eigenvalues of the A matrix (system poles).
///
/// A state matrix carrying a NaN or an infinity, and a solve that does not
/// converge, both yield an all-NaN array rather than a partially populated one:
/// the returned spectrum is either every eigenvalue or explicitly none.
///
/// @cite kailath1980 -- Kailath, "Linear Systems", 1980, Sec. 2.4
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
std::array<std::complex<Scalar>, NX> poles(const continuous_state_space<Scalar, NX, NU, NY>& sys)
{
    return detail::spectrum<Scalar, NX>(sys.A);
}

/// @copydoc poles
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
std::array<std::complex<Scalar>, NX> poles(const discrete_state_space<Scalar, NX, NU, NY>& sys)
{
    return detail::spectrum<Scalar, NX>(sys.A);
}

/// @brief Continuous system is provably stable iff all poles have negative real part.
///
/// A non-finite state matrix, and a spectrum the solver could not compute,
/// prove nothing and therefore answer false.
///
/// @cite franklin2015 -- Franklin et al., "Feedback Control of Dynamic Systems", 2015, Ch. 3 (Hurwitz criterion)
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

/// @brief Discrete system is provably stable iff all poles have magnitude < 1.
///
/// A non-finite state matrix, and a spectrum the solver could not compute,
/// prove nothing and therefore answer false.
///
/// @cite astrom1997 -- Astrom & Wittenmark, "Computer-Controlled Systems", 3rd ed., 1997, Sec. 3.3 (unit-circle criterion)
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

/// @brief Checks rank of controllability matrix [B, AB, A^2 B, ..., A^{n-1} B].
///
/// Returns true if rank equals NX (full state controllability). A non-finite
/// entry in either matrix answers false: the rank of a matrix containing a NaN
/// is not defined, so controllability is not provable from such a pair.
///
/// @cite kailath1980 -- Kailath, "Linear Systems", 1980, Sec. 2.4 (Kalman rank test)
template <typename Scalar, std::size_t NX, std::size_t NU>
bool is_controllable(const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NX, NU>& B)
{
    if(!A.allFinite() || !B.allFinite())
        return false;

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

/// @brief Checks rank of observability matrix [C; CA; CA^2; ...; CA^{n-1}].
///
/// Returns true if rank equals NX (full state observability). A non-finite
/// entry in either matrix answers false, for the same reason as the
/// controllability test: the rank of such a matrix is not defined.
///
/// @cite kailath1980 -- Kailath, "Linear Systems", 1980, Sec. 3.4 (dual Kalman rank test)
template <typename Scalar, std::size_t NX, std::size_t NY>
bool is_observable(const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NY, NX>& C)
{
    if(!A.allFinite() || !C.allFinite())
        return false;

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
//
// Both the three operands and the constructed closed-loop matrix are checked
// for finiteness: a finite operand triple can still overflow through the
// product B*K, and an infinity minus an infinity leaves a NaN in Acl, from
// which stability is not provable. Either condition answers false.
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
    for(int i = 0; i < n; ++i)
        if(std::abs(evals(i)) >= Scalar{1})
            return false;
    return true;
}

// is_stable_observer: checks if all eigenvalues of (A - L*C) are inside the unit circle.
// For discrete-time observer stability verification.
//
// Operands and the constructed observer matrix are checked for finiteness, on
// the same argument as the closed-loop test.
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
    for(int i = 0; i < n; ++i)
        if(std::abs(evals(i)) >= Scalar{1})
            return false;
    return true;
}

}

#endif
