#ifndef HPP_GUARD_CTRLPP_TESTS_FUZZ_CARE_QUAD_REFERENCE_H
#define HPP_GUARD_CTRLPP_TESTS_FUZZ_CARE_QUAD_REFERENCE_H

/// @brief An independent continuous Riccati reference at binary128, for two
/// states and one input.
///
/// INDEPENDENT IS THE WHOLE POINT, SO IT IS SPELLED OUT.
///
/// A test oracle that calls the library's own verdict function cannot fail when
/// that verdict is systematically wrong. This reference shares nothing with the
/// solver it checks:
///
///  * Different algorithm. Newton-Kleinman policy iteration. It never forms the
///    Hamiltonian matrix, never computes a Schur decomposition, never reorders
///    an invariant subspace and never extracts a solution from a basis. Each
///    step solves a continuous Lyapunov equation by dense elimination and
///    updates the gain. It additionally never inverts `A`, which matters here
///    and not in the discrete case: the continuous equation is well posed for a
///    singular `A`, so a reference that inverted it would refuse poses the
///    solver is required to answer.
///  * Different precision. `__float128`, a 113-bit significand against
///    binary64's 53, so the reference's own rounding is roughly `1.9e-34`
///    against the `2.2e-16` of the answer being judged -- eighteen orders below
///    anything being measured.
///  * Different arithmetic library: none at all. Every operation here is
///    written out over raw arrays, so no expression-template evaluation order
///    and no vectorization is shared with the solve.
///  * No transcendental calls, so no quad-precision math library is linked.
///    Only `+ - * /` and comparison, which lower to compiler runtime routines
///    both GCC and Clang already link. Convergence is tested on SQUARED norms
///    for the same reason -- a square root at this precision would be a library
///    call -- and so are the post-loop guards, which are sign tests on products
///    rather than comparisons of magnitudes.
///
/// Neither of the two libraries those paragraphs refuse is named anywhere in
/// this file, so searching it for their names reports what it uses rather than
/// what it talks about.
///
/// The iteration is started from the gain the returned solution implies, which
/// is stabilizing whenever the solver reports success, so it converges to the
/// nearby exact solution of the SAME pose. That is what makes the distance
/// between the two a forward error rather than a residual.
///
/// `quad`, `quad_reference_result`, `detail::quad_abs`, `detail::solve_4x4` and
/// `detail::quad_unit_roundoff` come from the discrete reference by including
/// it, and are not restated here. Two hand-copies of one shared numerical
/// primitive in this tree have already drifted apart, which is why primitives
/// here are called rather than spelled again. Lifting the shared pieces out
/// into a third header was considered and dropped: it is churn on a committed
/// target that buys no behavior.
///
/// @cite kleinman1968 -- Kleinman, "On an iterative technique for Riccati equation computations", 1968
/// @cite laub1979     -- Laub, "A Schur method for solving algebraic Riccati equations", 1979. Cited for the equation and for the characterization of its stabilizing solution; its METHOD is precisely what this file avoids.

#include "dare_quad_reference.h"

namespace ctrlpp::fuzz
{

/// @brief Refine a returned continuous Riccati solution at binary128.
///
/// @param A,B,Q,R the pose, in binary64 as the solver received it, widened
///        exactly (every binary64 value is a binary128 value, so widening the
///        pose introduces no error at all -- the reference solves the same
///        problem, not a nearby one).
/// @param P_seed the returned solution, used only to start the iteration.
/// @param max_steps the step budget. Reaching it without converging returns
///        `converged == false`.
///
/// Convergence is declared on the relative Frobenius step, tested squared, by
/// the same three clauses the discrete twin uses. The second is again the one
/// that actually stops the iteration: Newton-Kleinman converges quadratically
/// only down to the level where the elimination's own rounding dominates, and
/// then limit-cycles well above a naive multiple of unit roundoff. The
/// continuous iteration reaches that same precision floor for the same reason,
/// so a step that has stopped decreasing has converged as far as this precision
/// allows, and calling that failure would report failure on poses solved to
/// thirty digits.
///
/// Convergence alone is NOT sufficient here, and that is the one place this is
/// more than a transliteration of the discrete twin. The Stein operator is
/// singular exactly when some pair of closed-loop eigenvalues multiplies to one;
/// the continuous Lyapunov operator is singular whenever some pair SUMS to zero,
/// which every purely imaginary conjugate pair meets and which any marginally
/// stable closed loop sits arbitrarily close to. `detail::solve_4x4` detects only
/// an exactly zero pivot, so a near-singular operator returns a large meaningless
/// iterate whose step then stops decreasing -- and the second convergence clause
/// would report that as converged. The post-loop guards below are what turn it
/// into an abstention.
inline auto quad_refine_care(const double A[2][2], const double B[2], const double Q[2][2], double R, const double P_seed[2][2], int max_steps) -> quad_reference_result
{
    // 1024 times binary128's unit roundoff, carried over from the discrete twin
    // unchanged. The multiple is the accumulated rounding of the elimination and
    // the products at this size, not a fitted constant: the second convergence
    // clause is what actually stops the iteration on well-conditioned poses, and
    // this one only bounds how far it is allowed to keep trying.
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

    // The gain the seed implies: K = R^-1 B' X, one input so the inverted
    // quantity is a scalar -- and, unlike the discrete twin's R + B'XB, it is an
    // input of the problem rather than something the iterate can drive to zero.
    // The guard is kept for symmetry with that twin and to leave this header
    // self-contained; on the domain the consuming target admits it is
    // unreachable, because that target floors R strictly above zero.
    auto gain_from = [&](const quad S[2][2], quad K[2]) -> bool
    {
        if(Rq == quad{0})
            return false;
        K[0] = (Bq[0] * S[0][0] + Bq[1] * S[1][0]) / Rq;
        K[1] = (Bq[0] * S[0][1] + Bq[1] * S[1][1]) / Rq;
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

        // (I (x) Acl' + Acl' (x) I) vec(X) = vec(-W), column-major vec. The two
        // factors are the same two the discrete twin's Stein assembly
        // multiplies; here they are added under two Kronecker deltas instead.
        quad M[4][4];
        quad rhs[4];
        quad solution[4];
        for(int j2 = 0; j2 < 2; ++j2)
        {
            for(int i2 = 0; i2 < 2; ++i2)
            {
                const int row = i2 + 2 * j2;
                // THE SIGN IS THE POINT. The continuous Kleinman step is
                // Acl' X + X Acl + (Q + K' R K) = 0, so the right-hand side is
                // NEGATED relative to the discrete twin's rhs[row] = W[i2][j2].
                // Dropping the negation does not announce itself: it produces a
                // negative-definite X that converges just as cleanly, so the
                // convergence test cannot catch it. The post-loop definiteness
                // guard is what catches it.
                rhs[row] = -W[i2][j2];
                for(int j1 = 0; j1 < 2; ++j1)
                {
                    for(int i1 = 0; i1 < 2; ++i1)
                    {
                        const int col       = i1 + 2 * j1;
                        const quad left     = (j1 == j2) ? Acl[i1][i2] : quad{0};
                        const quad right    = (i1 == i2) ? Acl[j1][j2] : quad{0};
                        M[row][col]         = left + right;
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

    if(!result.converged)
        return result;

    // Three post-loop guards. Each can only WITHDRAW the verdict: a reference
    // that fails one reports itself as non-converged, which is an absence of
    // evidence and not a refutation of anything the solver said.

    // (1) Symmetry. The exact solution of a Lyapunov equation with a symmetric
    // right-hand side is symmetric, but the dense elimination is a general
    // solver and does not enforce it, so what survives is the elimination's own
    // accumulated rounding -- the same quantity the convergence tolerance above
    // already names. It is reused rather than a second bound being introduced,
    // and it is applied squared and relative for the same reason the step is.
    const quad asymmetry   = X[0][1] - X[1][0];
    const quad asymmetry_2 = asymmetry * asymmetry;
    quad magnitude_2       = quad{0};
    for(int i = 0; i < 2; ++i)
        for(int j = 0; j < 2; ++j)
            magnitude_2 += X[i][j] * X[i][j];
    if(asymmetry_2 > step_tolerance_2 * magnitude_2)
    {
        result.converged = false;
        return result;
    }

    // (2) Positive semi-definiteness of the reference's own answer. A real
    // symmetric matrix is positive semi-definite exactly when ALL its principal
    // minors are non-negative, which at this dimension is the two diagonal
    // entries and the determinant -- the leading minors alone would not do, as
    // diag(0, -1) passes those and is indefinite. The test is a comparison of
    // products against zero: exact, with no square root and no tolerance of its
    // own. This is the guard that catches a dropped negation on the right-hand
    // side above, which the convergence test cannot see.
    const quad solution_determinant = X[0][0] * X[1][1] - X[0][1] * X[1][0];
    if(X[0][0] < quad{0} || X[1][1] < quad{0} || solution_determinant < quad{0})
    {
        result.converged = false;
        return result;
    }

    // (3) Strict stability of the closed loop the FINAL solution implies. This
    // is what separates the stabilizing solution from the other real symmetric
    // solutions of the same equation, and it is also what catches the
    // near-singular Lyapunov operator described above: the large meaningless
    // iterate it returns does not carry a stable closed loop. For a real 2-by-2
    // the characteristic polynomial is lambda^2 - trace*lambda + determinant, so
    // both roots lie strictly left of the imaginary axis exactly when the trace
    // is strictly negative and the determinant strictly positive. That is the
    // Routh-Hurwitz condition at this dimension: two sign tests, no eigensolver,
    // no square root and no tolerance.
    if(!gain_from(X, K))
    {
        result.converged = false;
        return result;
    }
    quad Acl_final[2][2];
    for(int i = 0; i < 2; ++i)
        for(int j = 0; j < 2; ++j)
            Acl_final[i][j] = Aq[i][j] - Bq[i] * K[j];

    const quad closed_loop_trace       = Acl_final[0][0] + Acl_final[1][1];
    const quad closed_loop_determinant = Acl_final[0][0] * Acl_final[1][1] - Acl_final[0][1] * Acl_final[1][0];
    if(!(closed_loop_trace < quad{0}) || !(closed_loop_determinant > quad{0}))
        result.converged = false;

    return result;
}

}

#endif
