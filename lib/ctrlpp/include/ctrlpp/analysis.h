#ifndef HPP_GUARD_CTRLPP_ANALYSIS_H
#define HPP_GUARD_CTRLPP_ANALYSIS_H

#include "ctrlpp/state_space.h"

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>

namespace ctrlpp {

namespace detail {

template<typename Scalar, std::size_t N, typename Policy>
concept EigenvaluePolicy = LinalgPolicy<Policy> && requires(
    typename Policy::template matrix_type<Scalar, N, N> A
) {
    { Policy::eigenvalues(A) } -> std::same_as<std::array<std::complex<Scalar>, N>>;
};

}

// poles: returns eigenvalues of the A matrix (system poles).
template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, LinalgPolicy Policy>
    requires detail::EigenvaluePolicy<Scalar, NX, Policy>
auto poles(const ContinuousStateSpace<Scalar, NX, NU, NY, Policy>& sys)
    -> std::array<std::complex<Scalar>, NX>
{
    return Policy::eigenvalues(sys.A);
}

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, LinalgPolicy Policy>
    requires detail::EigenvaluePolicy<Scalar, NX, Policy>
auto poles(const DiscreteStateSpace<Scalar, NX, NU, NY, Policy>& sys)
    -> std::array<std::complex<Scalar>, NX>
{
    return Policy::eigenvalues(sys.A);
}

// is_stable: continuous system is stable iff all poles have negative real part.
template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, LinalgPolicy Policy>
    requires detail::EigenvaluePolicy<Scalar, NX, Policy>
auto is_stable(const ContinuousStateSpace<Scalar, NX, NU, NY, Policy>& sys) -> bool
{
    auto p = poles(sys);
    for (const auto& pole : p)
        if (pole.real() >= Scalar{0})
            return false;
    return true;
}

// is_stable: discrete system is stable iff all poles have magnitude < 1.
template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, LinalgPolicy Policy>
    requires detail::EigenvaluePolicy<Scalar, NX, Policy>
auto is_stable(const DiscreteStateSpace<Scalar, NX, NU, NY, Policy>& sys) -> bool
{
    auto p = poles(sys);
    for (const auto& pole : p)
        if (std::abs(pole) >= Scalar{1})
            return false;
    return true;
}

// is_controllable: checks rank of controllability matrix [B, AB, A^2 B, ..., A^{n-1} B].
// Returns true if rank equals NX (full state controllability).
// Requires Policy to provide rank() for an NX x (NX*NU) matrix.
template<typename Scalar, std::size_t NX, std::size_t NU, LinalgPolicy Policy>
    requires requires(typename Policy::template matrix_type<Scalar, NX, NX * NU> M) {
        { Policy::rank(M) } -> std::convertible_to<std::size_t>;
    }
auto is_controllable(
    const typename Policy::template matrix_type<Scalar, NX, NX>& A,
    const typename Policy::template matrix_type<Scalar, NX, NU>& B) -> bool
{
    using MatB = typename Policy::template matrix_type<Scalar, NX, NU>;
    using MatC = typename Policy::template matrix_type<Scalar, NX, NX * NU>;

    MatC C{};
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);

    // Set first block: C[:, 0:NU] = B
    for (int r = 0; r < nx; ++r)
        for (int c = 0; c < nu; ++c)
            Policy::set_element(C, static_cast<std::size_t>(r),
                                static_cast<std::size_t>(c),
                                Policy::get_element(B, static_cast<std::size_t>(r),
                                                    static_cast<std::size_t>(c)));

    MatB AkB = B;
    for (std::size_t k = 1; k < NX; ++k) {
        AkB = Policy::multiply(A, AkB);
        for (int r = 0; r < nx; ++r)
            for (int c = 0; c < nu; ++c)
                Policy::set_element(C, static_cast<std::size_t>(r),
                                    k * NU + static_cast<std::size_t>(c),
                                    Policy::get_element(AkB, static_cast<std::size_t>(r),
                                                        static_cast<std::size_t>(c)));
    }

    return Policy::rank(C) == NX;
}

// is_observable: checks rank of observability matrix [C; CA; CA^2; ...; CA^{n-1}].
// Returns true if rank equals NX (full state observability).
template<typename Scalar, std::size_t NX, std::size_t NY, LinalgPolicy Policy>
    requires requires(typename Policy::template matrix_type<Scalar, NX * NY, NX> M) {
        { Policy::rank(M) } -> std::convertible_to<std::size_t>;
    }
auto is_observable(
    const typename Policy::template matrix_type<Scalar, NX, NX>& A,
    const typename Policy::template matrix_type<Scalar, NY, NX>& C) -> bool
{
    using MatC = typename Policy::template matrix_type<Scalar, NY, NX>;
    using MatO = typename Policy::template matrix_type<Scalar, NX * NY, NX>;

    MatO O{};
    constexpr int nx = static_cast<int>(NX);
    constexpr int ny = static_cast<int>(NY);

    // Set first block: O[0:NY, :] = C
    for (int r = 0; r < ny; ++r)
        for (int c = 0; c < nx; ++c)
            Policy::set_element(O, static_cast<std::size_t>(r),
                                static_cast<std::size_t>(c),
                                Policy::get_element(C, static_cast<std::size_t>(r),
                                                    static_cast<std::size_t>(c)));

    MatC CAk = C;
    for (std::size_t k = 1; k < NX; ++k) {
        CAk = Policy::multiply(CAk, A);
        for (int r = 0; r < ny; ++r)
            for (int c = 0; c < nx; ++c)
                Policy::set_element(O, k * NY + static_cast<std::size_t>(r),
                                    static_cast<std::size_t>(c),
                                    Policy::get_element(CAk, static_cast<std::size_t>(r),
                                                        static_cast<std::size_t>(c)));
    }

    return Policy::rank(O) == NX;
}

// is_stable_closed_loop: checks if all eigenvalues of (A - B*K) are inside the unit circle.
// For discrete-time closed-loop stability verification.
template<typename Scalar, std::size_t NX, std::size_t NU, LinalgPolicy Policy>
    requires detail::EigenvaluePolicy<Scalar, NX, Policy>
auto is_stable_closed_loop(
    const typename Policy::template matrix_type<Scalar, NX, NX>& A,
    const typename Policy::template matrix_type<Scalar, NX, NU>& B,
    const typename Policy::template matrix_type<Scalar, NU, NX>& K) -> bool
{
    auto Acl = Policy::subtract(A, Policy::multiply(B, K));
    auto evals = Policy::eigenvalues(Acl);
    for (const auto& ev : evals)
        if (std::abs(ev) >= Scalar{1})
            return false;
    return true;
}

// is_stable_observer: checks if all eigenvalues of (A - L*C) are inside the unit circle.
// For discrete-time observer stability verification.
template<typename Scalar, std::size_t NX, std::size_t NY, LinalgPolicy Policy>
    requires detail::EigenvaluePolicy<Scalar, NX, Policy>
auto is_stable_observer(
    const typename Policy::template matrix_type<Scalar, NX, NX>& A,
    const typename Policy::template matrix_type<Scalar, NX, NY>& L,
    const typename Policy::template matrix_type<Scalar, NY, NX>& C) -> bool
{
    auto Aobs = Policy::subtract(A, Policy::multiply(L, C));
    auto evals = Policy::eigenvalues(Aobs);
    for (const auto& ev : evals)
        if (std::abs(ev) >= Scalar{1})
            return false;
    return true;
}

}

#endif
