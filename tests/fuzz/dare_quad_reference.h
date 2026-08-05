#ifndef HPP_GUARD_CTRLPP_TESTS_FUZZ_DARE_QUAD_REFERENCE_H
#define HPP_GUARD_CTRLPP_TESTS_FUZZ_DARE_QUAD_REFERENCE_H

/// @brief An independent discrete Riccati reference at binary128, for two states
/// and one input.
///
/// INDEPENDENT IS THE WHOLE POINT, SO IT IS SPELLED OUT.
///
/// A test oracle that calls the library's own verdict function cannot fail when
/// that verdict is systematically wrong. This reference shares nothing with the
/// solver it checks:
///
///  * Different algorithm. Newton-Kleinman policy iteration. It never forms the
///    symplectic matrix, never computes a Schur decomposition, never reorders an
///    invariant subspace and never extracts a solution from a basis. Each step
///    solves a Stein equation by dense elimination and updates the gain.
///  * Different precision. `__float128`, a 113-bit significand against
///    binary64's 53, so the reference's own rounding is roughly `1.9e-34`
///    against the `2.2e-16` of the answer being judged -- eighteen orders below
///    anything being measured.
///  * Different arithmetic library. No Eigen. Every operation here is written
///    out, so no expression-template evaluation order and no vectorization is
///    shared with the solve.
///  * No transcendental calls, so no `libquadmath`. Only `+ - * /` and
///    comparison, which lower to compiler runtime routines both GCC and Clang
///    already link. Convergence is tested on SQUARED norms for the same reason:
///    a square root at this precision would be a library call.
///
/// The iteration is started from the gain the returned solution implies, which
/// is stabilizing whenever the solver reports success, so it converges to the
/// nearby exact solution of the SAME pose. That is what makes the distance
/// between the two a forward error rather than a residual.
///
/// @cite kleinman1968 -- Kleinman, "On an iterative technique for Riccati equation computations", 1968
/// @cite hewer1971    -- Hewer, "An iterative technique for the computation of the steady state gains for the discrete optimal regulator", 1971

#include <cstddef>

namespace ctrlpp::fuzz
{

using quad = __float128;

/// @brief What a reference run produced, and whether it is usable as one.
///
/// `converged` false is an ABSENCE OF EVIDENCE, not a refutation. A reference
/// that has not converged is not a reference, and a caller must not read a
/// distance to it as an error.
struct quad_reference_result
{
    quad P[2][2]{};
    bool converged{};
    int steps{};
};

namespace detail
{

inline auto quad_abs(quad x) -> quad
{
    return x < quad{0} ? -x : x;
}

/// @brief binary128's unit roundoff, 2^-112.
///
/// Built by halving, so the header needs neither a quad literal suffix nor a
/// quad-precision math header -- either would make this reference depend on a
/// library it exists to be independent of.
inline auto quad_unit_roundoff() -> quad
{
    quad e = quad{1};
    for(int i = 0; i < 112; ++i)
        e /= quad{2};
    return e;
}

/// @brief Solve a 4x4 system by Gaussian elimination with partial pivoting.
///
/// Returns false when a pivot is exactly zero, which is the only singularity
/// this routine can detect without a scale to compare against; the caller treats
/// that as non-convergence rather than as a result.
inline auto solve_4x4(quad M[4][4], quad rhs[4], quad out[4]) -> bool
{
    for(int col = 0; col < 4; ++col)
    {
        int pivot = col;
        for(int row = col + 1; row < 4; ++row)
        {
            if(quad_abs(M[row][col]) > quad_abs(M[pivot][col]))
                pivot = row;
        }
        if(M[pivot][col] == quad{0})
            return false;
        if(pivot != col)
        {
            for(int k = 0; k < 4; ++k)
            {
                const quad swap = M[col][k];
                M[col][k]       = M[pivot][k];
                M[pivot][k]     = swap;
            }
            const quad swap = rhs[col];
            rhs[col]        = rhs[pivot];
            rhs[pivot]      = swap;
        }
        for(int row = col + 1; row < 4; ++row)
        {
            const quad factor = M[row][col] / M[col][col];
            if(factor == quad{0})
                continue;
            for(int k = col; k < 4; ++k)
                M[row][k] -= factor * M[col][k];
            rhs[row] -= factor * rhs[col];
        }
    }
    for(int row = 3; row >= 0; --row)
    {
        quad accumulated = rhs[row];
        for(int k = row + 1; k < 4; ++k)
            accumulated -= M[row][k] * out[k];
        out[row] = accumulated / M[row][row];
    }
    return true;
}

}

/// @brief Refine a returned discrete Riccati solution at binary128.
///
/// @param A,B,Q,R the pose, in binary64 as the solver received it, widened
///        exactly (every binary64 value is a binary128 value, so widening the
///        pose introduces no error at all -- the reference solves the same
///        problem, not a nearby one).
/// @param P_seed the returned solution, used only to start the iteration.
/// @param max_steps the step budget. Reaching it without converging returns
///        `converged == false`.
///
/// Convergence is declared on the relative Frobenius step, tested squared. Two
/// clauses, and the second is load-bearing: Newton-Kleinman converges
/// quadratically only down to the level where the elimination's own rounding
/// dominates, and then limit-cycles well above a naive multiple of unit
/// roundoff. A step that has stopped decreasing has converged as far as this
/// precision allows, and treating that as failure would report failure on poses
/// solved to thirty digits.
inline auto quad_refine_dare(const double A[2][2], const double B[2], const double Q[2][2], double R, const double P_seed[2][2], int max_steps) -> quad_reference_result
{
    // 1024 times binary128's unit roundoff. The multiple is the accumulated
    // rounding of the elimination and the products at this size, not a fitted
    // constant: the second convergence clause is what actually stops the
    // iteration on well-conditioned poses, and this one only bounds how far it
    // is allowed to keep trying.
    const quad quad_epsilon     = detail::quad_unit_roundoff();
    const quad step_tolerance   = quad{1024} * quad_epsilon;
    const quad step_tolerance_2 = step_tolerance * step_tolerance;

    quad Aq[2][2];
    quad Bq[2];
    quad Qq[2][2];
    for(int i = 0; i < 2; ++i)
    {
        Bq[i] = quad{B[i]};
        for(int j = 0; j < 2; ++j)
        {
            Aq[i][j] = quad{A[i][j]};
            Qq[i][j] = quad{Q[i][j]};
        }
    }
    const quad Rq = quad{R};

    quad X[2][2];
    for(int i = 0; i < 2; ++i)
        for(int j = 0; j < 2; ++j)
            X[i][j] = quad{P_seed[i][j]};

    // The gain the seed implies: K = (R + B'XB)^-1 B'XA, one input so the
    // inverted quantity is a scalar.
    auto gain_from = [&](const quad S[2][2], quad K[2]) -> bool
    {
        quad BtS[2] = {Bq[0] * S[0][0] + Bq[1] * S[1][0], Bq[0] * S[0][1] + Bq[1] * S[1][1]};
        const quad denominator = Rq + BtS[0] * Bq[0] + BtS[1] * Bq[1];
        if(denominator == quad{0})
            return false;
        K[0] = (BtS[0] * Aq[0][0] + BtS[1] * Aq[1][0]) / denominator;
        K[1] = (BtS[0] * Aq[0][1] + BtS[1] * Aq[1][1]) / denominator;
        return true;
    };

    quad K[2];
    if(!gain_from(X, K))
        return {};

    quad_reference_result result;
    quad previous_step_2 = quad{-1};

    for(int step = 1; step <= max_steps; ++step)
    {
        // A_cl = A - B K
        quad Acl[2][2];
        for(int i = 0; i < 2; ++i)
            for(int j = 0; j < 2; ++j)
                Acl[i][j] = Aq[i][j] - Bq[i] * K[j];

        // W = Q + K' R K
        quad W[2][2];
        for(int i = 0; i < 2; ++i)
            for(int j = 0; j < 2; ++j)
                W[i][j] = Qq[i][j] + K[i] * Rq * K[j];

        // (I - Acl' (x) Acl') vec(X) = vec(W), column-major vec.
        quad M[4][4];
        quad rhs[4];
        quad solution[4];
        for(int j2 = 0; j2 < 2; ++j2)
        {
            for(int i2 = 0; i2 < 2; ++i2)
            {
                const int row = i2 + 2 * j2;
                rhs[row]      = W[i2][j2];
                for(int j1 = 0; j1 < 2; ++j1)
                {
                    for(int i1 = 0; i1 < 2; ++i1)
                    {
                        const int col = i1 + 2 * j1;
                        const quad identity = (i1 == i2 && j1 == j2) ? quad{1} : quad{0};
                        M[row][col]         = identity - Acl[j1][j2] * Acl[i1][i2];
                    }
                }
            }
        }
        if(!detail::solve_4x4(M, rhs, solution))
            return result;

        quad X_next[2][2];
        X_next[0][0] = solution[0];
        X_next[1][0] = solution[1];
        X_next[0][1] = solution[2];
        X_next[1][1] = solution[3];

        quad difference_2 = quad{0};
        quad magnitude_2  = quad{0};
        for(int i = 0; i < 2; ++i)
        {
            for(int j = 0; j < 2; ++j)
            {
                const quad d = X_next[i][j] - X[i][j];
                difference_2 += d * d;
                magnitude_2 += X_next[i][j] * X_next[i][j];
            }
        }

        for(int i = 0; i < 2; ++i)
            for(int j = 0; j < 2; ++j)
                X[i][j] = X_next[i][j];

        result.steps = step;

        if(magnitude_2 == quad{0})
        {
            result.converged = (difference_2 == quad{0});
            break;
        }

        const quad step_2 = difference_2 / magnitude_2;
        if(step_2 <= step_tolerance_2)
        {
            result.converged = true;
            break;
        }
        // The step stopped decreasing: the iteration has reached the floor this
        // precision supports and further steps only move rounding around.
        if(previous_step_2 >= quad{0} && step_2 >= previous_step_2)
        {
            result.converged = true;
            break;
        }
        previous_step_2 = step_2;

        if(!gain_from(X, K))
            return result;
    }

    for(int i = 0; i < 2; ++i)
        for(int j = 0; j < 2; ++j)
            result.P[i][j] = X[i][j];
    return result;
}

/// @brief Relative Frobenius distance from a returned solution to the reference,
/// SQUARED, so the comparison needs no square root at binary128.
///
/// Returns false when the reference's own magnitude is zero, which carries no
/// scale to be relative to.
inline auto quad_relative_distance_squared(const double P[2][2], const quad P_reference[2][2], quad &out) -> bool
{
    quad difference_2 = quad{0};
    quad magnitude_2  = quad{0};
    for(int i = 0; i < 2; ++i)
    {
        for(int j = 0; j < 2; ++j)
        {
            const quad d = quad{P[i][j]} - P_reference[i][j];
            difference_2 += d * d;
            magnitude_2 += P_reference[i][j] * P_reference[i][j];
        }
    }
    if(magnitude_2 == quad{0})
        return false;
    out = difference_2 / magnitude_2;
    return true;
}

}

#endif
