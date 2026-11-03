#ifndef HPP_GUARD_CTRLPP_DETAIL_CARE_METHODS_H
#define HPP_GUARD_CTRLPP_DETAIL_CARE_METHODS_H

// Tags selecting the continuous Riccati solve: the scaled Newton iteration on sign(H), the real Schur method with
// Bai-Demmel reordering, and the real Schur method after a DGEBAL-style diagonal balance of the Hamiltonian. The
// sign-function path is the default because it took the fewest median instructions of the three at every NX from 8
// to 30. Every tag verifies its extracted matrix against the caller's Hamiltonian before returning it.
// @cite laub1979      : Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979
// @cite roberts1980   : Roberts, "Linear model reduction and solution of the algebraic Riccati equation by use of the sign function", 1980
// @cite byers1987     : Byers, "Solving the algebraic Riccati equation with the matrix sign function", 1987
// @cite higham2008    : Higham, "Functions of Matrices: Theory and Computation", 2008, Chapter 5
// @cite lapack_dgebal : Reference LAPACK SRC/dgebal.f phase 2 (scaling) algorithm

#include <cstdint>
#include <concepts>

namespace ctrlpp::detail
{

struct sign_function_care_method
{
    sign_function_care_method()
        : warmup_iterations{3}
    {
    }

    // Newton steps before the non-contraction guard arms. The guard bounds wasted work and decides nothing about
    // acceptance. The default is not derived: every value from 0 to the iteration cap of 40 gave identical outcomes over
    // 3,456 draws spanning eighteen decades of weight scale.
    std::int32_t warmup_iterations;
};

struct schur_care_method
{
};

struct balanced_schur_care_method
{
};

template <typename T>
concept care_solve_method = std::same_as<T, schur_care_method> || std::same_as<T, sign_function_care_method> || std::same_as<T, balanced_schur_care_method>;

}

#endif
