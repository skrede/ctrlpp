#ifndef HPP_GUARD_CTRLPP_DETAIL_SCHUR_REORDER_H
#define HPP_GUARD_CTRLPP_DETAIL_SCHUR_REORDER_H

/// @brief Real-Schur block reordering via Bai-Demmel 1993 predicated swap kernel.
///
/// The primitive supports all four adjacent-block swap configurations
/// (1x1/1x1, 1x1/2x2, 2x2/1x1, 2x2/2x2) on a real quasi-triangular T with
/// orthogonal basis U produced by Eigen::RealSchur. A caller-supplied
/// predicate selects which eigenvalues to move to the leading position;
/// DLANV2 Murnaghan standardisation restores the canonical form of every
/// touched 2x2 block. Every threshold is derived from
/// std::numeric_limits<Scalar>::epsilon() * matrix_scale, matching the
/// LAPACK DTREXC / DLASY2 contract in pure real arithmetic on top of
/// fixed-size Eigen matrices.
///
/// @cite bai_demmel_1993 -- Bai & Demmel, "On swapping diagonal blocks in real Schur form", 1993
/// @cite laub1979       -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979

#include "ctrlpp/detail/quasi_triangular.h"

#include <Eigen/QR>
#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <complex>
#include <utility>
#include <concepts>
#include <algorithm>
#include <type_traits>

namespace ctrlpp::detail
{

// --- Swap conditioning policy tags ---
//
// Selects the test used to reject an ill-conditioned Bai-Demmel swap.
//  * pivot_ratio_conditioning: |R(last,last)| / |R(0,0)| from the rank-revealing
//    ColPivHouseholderQR on the Sylvester operator (default, zero extra cost).
//  * hager_higham_conditioning: Hager-Higham 1-norm estimate of K^-1 (tighter,
//    3-5 extra 4x4 matvecs per swap). Not yet implemented -- compile-error stub.

struct pivot_ratio_conditioning
{
};

struct hager_higham_conditioning
{
};

/// Concept satisfied by the two conditioning policy tag types. Public API
/// templates (dare, care) constrain their trailing Cond parameter on this so
/// that overload resolution rejects foreign types in that slot -- preventing
/// ambiguity with overloads that accept an additional matrix argument.
template <typename T>
concept conditioning_policy =
    std::same_as<T, pivot_ratio_conditioning>
 || std::same_as<T, hager_higham_conditioning>;

// --- Precision-templated constexpr multipliers (no bare literals) ---
//
// detail_constants wraps LAPACK-derived multipliers behind named constexpr
// functions so that the bodies of the reorder helpers contain no numeric
// literals other than structural identity elements.

namespace detail_constants
{

// DLANV2 real-vs-complex-pair decision multiplier (matches LAPACK DLANV2 'MULTPL' = 4).
template <typename Scalar>
constexpr auto dlanv2_multpl() noexcept -> Scalar { return Scalar{4}; }

// LAPACK SLAEXC post-swap residual threshold multiplier (default 10).
template <typename Scalar>
constexpr auto rejection_multiplier() noexcept -> Scalar { return Scalar{10}; }

}

// --- Driver result POD ---
//
//  * placed               : number of predicate-matching eigenvalues moved to the top.
//  * complete             : false if any swap was rejected by the conditioning test.
//  * subspace_separation  : min pivot ratio across all accepted swaps (LAPACK SEP analogue).

template <typename Scalar>
struct reorder_result
{
    int    placed{0};
    bool   complete{true};
    Scalar subspace_separation{Scalar{1}};
};

// --- Forward declarations (definitions follow) ---

template <typename Scalar, int N>
auto standardize_2x2_block(Eigen::Matrix<Scalar, N, N>& T,
                           Eigen::Matrix<Scalar, N, N>& U,
                           int p) -> void;

template <typename Scalar, int N, typename Cond>
auto swap_real_schur_blocks(Eigen::Matrix<Scalar, N, N>& T,
                            Eigen::Matrix<Scalar, N, N>& U,
                            int p, int n1, int n2,
                            Cond /*tag*/,
                            Scalar& pivot_ratio_out) -> bool;

template <typename Scalar, int N, typename Predicate,
          typename Cond = pivot_ratio_conditioning>
auto reorder_real_schur(Eigen::Matrix<Scalar, N, N>& T,
                        Eigen::Matrix<Scalar, N, N>& U,
                        Predicate&&                   predicate,
                        Cond                          tag = {})
    -> reorder_result<Scalar>
{
    static_assert(std::is_same_v<Cond, pivot_ratio_conditioning>
                  || std::is_same_v<Cond, hager_higham_conditioning>,
                  "Cond must be pivot_ratio_conditioning or hager_higham_conditioning");
    if constexpr (std::is_same_v<Cond, hager_higham_conditioning>)
    {
        static_assert(!std::is_same_v<Cond, hager_higham_conditioning>,
                      "hager_higham_conditioning not yet implemented -- use pivot_ratio_conditioning");
    }

    reorder_result<Scalar> r{};
    r.complete = true;
    r.subspace_separation = Scalar{1};

    const Scalar scale_T = T.cwiseAbs().maxCoeff();

    // Block structure and block eigenvalues come from the shared
    // quasi-triangular primitives rather than from a copy specialized to this
    // caller: the significance test on the subdiagonal entry and the
    // trace-and-determinant eigenvalue form are the same two questions the
    // continuous acceptance rule asks of its own factor, and two copies of a
    // predicate are two contracts.
    auto block_size_at = [&](int pos) -> int
    {
        return quasi_triangular_block_size<Scalar, N>(T, pos, scale_T);
    };

    auto block_matches = [&](int pos, int n_block) -> bool
    {
        // `first` is the root with the larger real part on a real pair and the
        // shared real part on a conjugate pair, so one query answers both. For
        // a conjugate pair |lambda| and sign(Re lambda) are identical across
        // the two roots in any case.
        return predicate(
            quasi_triangular_block_spectrum<Scalar, N>(T, pos, n_block).first);
    };

    int placed_pos = 0;
    while (placed_pos < N)
    {
        // Scan right from placed_pos for the first predicate-matching block.
        int scan = placed_pos;
        int scan_nb = 1;
        bool found = false;
        while (scan < N)
        {
            const int nb = block_size_at(scan);
            if (block_matches(scan, nb))
            {
                scan_nb = nb;
                found = true;
                break;
            }
            scan += nb;
        }

        if (!found)
            break;

        // Bubble the found block from `scan` leftward to `placed_pos`.
        int cur = scan;
        while (cur > placed_pos)
        {
            // Determine the block immediately to the left of `cur`. A 2x2 block
            // ending at cur - 1 starts at cur - 2, so this is the same
            // significance question asked one block earlier.
            const int left_nb = (cur >= 2) ? block_size_at(cur - 2) : 1;
            const int left_pos = cur - left_nb;

            Scalar pivot_ratio = Scalar{1};
            const bool ok = swap_real_schur_blocks<Scalar, N, Cond>(
                T, U, left_pos, left_nb, scan_nb, tag, pivot_ratio);

            if (!ok)
            {
                r.complete = false;
                r.subspace_separation = std::min(r.subspace_separation, pivot_ratio);
                break;  // Abort this bubble; the result reports a partial reorder.
            }

            r.subspace_separation = std::min(r.subspace_separation, pivot_ratio);

            // Standardize any 2x2 block touched by the swap.
            if (scan_nb == 2)
                standardize_2x2_block<Scalar, N>(T, U, left_pos);
            if (left_nb == 2)
                standardize_2x2_block<Scalar, N>(T, U, left_pos + scan_nb);

            cur = left_pos;
            // Re-probe scan_nb: standardisation can split a 2x2 block into
            // two 1x1 blocks (complex pair degenerating to real eigenvalues).
            scan_nb = block_size_at(cur);
            if (!block_matches(cur, scan_nb))
            {
                // The bubbling block lost its predicate match due to a split;
                // abort; the outer scan will find the next candidate.
                break;
            }
        }

        if (cur == placed_pos)
        {
            r.placed += scan_nb;
            placed_pos += scan_nb;
        }
        else
        {
            // Bubble aborted before reaching placed_pos; resume scanning
            // past the stuck block (progress guarantees termination).
            placed_pos = cur + block_size_at(cur);
        }
    }

    return r;
}

// --- DLANV2 Murnaghan 2x2 block standardization ---
//
// Restores the canonical form of a 2x2 diagonal block at T(p:p+2, p:p+2):
//   * real eigenvalue pair  -> block becomes upper triangular (C = 0).
//   * complex-conjugate pair -> equal diagonals, opposite-sign off-diagonals,
//                               |b|*|c| = |lambda_im|^2.
//
// Accumulates the Givens rotation into U so that T = U S U^T holds with the
// incoming Schur basis S. Matches LAPACK DLANV2 branch structure (C=0,
// B=0, A=D, general) verbatim; the real-vs-complex-pair decision threshold
// MULTPL * eps is lifted into detail_constants::dlanv2_multpl<Scalar>().

template <typename Scalar, int N>
auto standardize_2x2_block(Eigen::Matrix<Scalar, N, N>& T,
                           Eigen::Matrix<Scalar, N, N>& U,
                           int p) -> void
{
    using std::abs;
    using std::sqrt;
    using std::hypot;
    using std::copysign;

    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    const Scalar multpl = detail_constants::dlanv2_multpl<Scalar>();

    const Scalar a = T(p,     p);
    const Scalar b = T(p,     p + 1);
    const Scalar c = T(p + 1, p);
    const Scalar d = T(p + 1, p + 1);

    Scalar cs = Scalar{1};
    Scalar sn = Scalar{0};
    bool   enforce_equal_diag = false;
    Scalar ad_mean = Scalar{0};

    // Branch 1: C = 0 -- block already upper triangular. No-op.
    if (c == Scalar{0})
    {
        cs = Scalar{1};
        sn = Scalar{0};
    }
    // Branch 2: B = 0 -- swap rows/cols to make strictly upper triangular.
    else if (b == Scalar{0})
    {
        cs = Scalar{0};
        sn = Scalar{1};
    }
    // Branch 3: A == D and sign(B) != sign(C) -- already canonical.
    else if (a == d && ((b > Scalar{0}) != (c > Scalar{0})))
    {
        cs = Scalar{1};
        sn = Scalar{0};
    }
    else
    {
        const Scalar temp  = a - d;
        const Scalar p_val = temp / Scalar{2};
        const Scalar bcmax = std::max(abs(b), abs(c));
        const Scalar sgn_b = (b < Scalar{0}) ? Scalar{-1} : Scalar{1};
        const Scalar sgn_c = (c < Scalar{0}) ? Scalar{-1} : Scalar{1};
        const Scalar bcmis = std::min(abs(b), abs(c)) * sgn_b * sgn_c;
        const Scalar scl   = std::max(abs(p_val), bcmax);
        const Scalar z_val = (p_val / scl) * p_val + (bcmax / scl) * bcmis;

        if (z_val >= multpl * eps)
        {
            // Real eigenvalues branch.
            Scalar z   = p_val + copysign(sqrt(scl) * sqrt(z_val), p_val);
            Scalar tau = hypot(c, z);
            cs = z / tau;
            sn = c / tau;
        }
        else
        {
            // Complex-conjugate pair branch.
            const Scalar sigma   = b + c;
            const Scalar sgn_sig = (sigma < Scalar{0}) ? Scalar{-1} : Scalar{1};
            Scalar tau = hypot(sigma, temp);
            cs = sqrt((Scalar{1} + abs(sigma) / tau) / Scalar{2});
            sn = -(p_val / (tau * cs)) * sgn_sig;
            enforce_equal_diag = true;
            ad_mean = (a + d) / Scalar{2};
        }
    }

    // Apply Givens similarity: left on rows (p, p+1) of T.
    for (int j = 0; j < N; ++j)
    {
        const Scalar t0 = T(p,     j);
        const Scalar t1 = T(p + 1, j);
        T(p,     j) =  cs * t0 + sn * t1;
        T(p + 1, j) = -sn * t0 + cs * t1;
    }
    // Right on cols (p, p+1) of T and U (U accumulates the right-rotation).
    for (int i = 0; i < N; ++i)
    {
        const Scalar t0 = T(i, p);
        const Scalar t1 = T(i, p + 1);
        T(i, p)     =  cs * t0 + sn * t1;
        T(i, p + 1) = -sn * t0 + cs * t1;

        const Scalar u0 = U(i, p);
        const Scalar u1 = U(i, p + 1);
        U(i, p)     =  cs * u0 + sn * u1;
        U(i, p + 1) = -sn * u0 + cs * u1;
    }

    if (enforce_equal_diag)
    {
        // DLANV2 complex-pair branch: force exactly equal diagonals to
        // suppress floating-point drift in the sigma-derived rotation.
        T(p,     p)     = ad_mean;
        T(p + 1, p + 1) = ad_mean;
    }
}

// --- 1x1/1x1 closed-form swap (always succeeds per Bai-Demmel 1993, p.79) ---
//
// Solves the scalar Sylvester equation t11 * x - x * t22 = t12 for x, then
// applies the overflow-safe Givens rotation (Golub-Van Loan Alg. 5.1.3) to
// T (two-sided) and U (right-multiply). When eigenvalues coincide within
// eps * max(|t11|, |t22|) the swap reduces to a no-op (LAPACK convention).

template <typename Scalar, int N>
auto swap_real_schur_1x1(Eigen::Matrix<Scalar, N, N>& T,
                         Eigen::Matrix<Scalar, N, N>& U,
                         int p, Scalar& pivot_ratio_out) -> bool
{
    using std::abs;
    using std::sqrt;
    using std::hypot;

    const Scalar eps = std::numeric_limits<Scalar>::epsilon();

    const Scalar t11 = T(p,     p);
    const Scalar t22 = T(p + 1, p + 1);
    const Scalar t12 = T(p,     p + 1);
    const Scalar diff = t11 - t22;

    // Coinciding eigenvalues -> no-op (LAPACK: swap is trivial).
    if (abs(diff) < eps * std::max(abs(t11), abs(t22)))
    {
        pivot_ratio_out = Scalar{1};
        return true;
    }

    // G = [[cs, -sn], [sn, cs]] chosen so that G^T * [[t11, t12], [0, t22]] * G
    // has zero at position (1, 0). Solving for s/c == -diff/t12 gives
    // (cs, sn) = (t12, -diff) / sqrt(t12^2 + diff^2); std::hypot avoids
    // overflow/underflow for ill-scaled inputs.
    const Scalar r = hypot(t12, diff);
    const Scalar cs = t12 / r;
    const Scalar sn = -diff / r;

    // Similarity on rows (p, p+1) of T.
    for (int j = 0; j < N; ++j)
    {
        const Scalar a0 = T(p,     j);
        const Scalar a1 = T(p + 1, j);
        T(p,     j) =  cs * a0 + sn * a1;
        T(p + 1, j) = -sn * a0 + cs * a1;
    }
    // Similarity on cols (p, p+1) of T; U accumulates right-rotations.
    for (int i = 0; i < N; ++i)
    {
        const Scalar a0 = T(i, p);
        const Scalar a1 = T(i, p + 1);
        T(i, p)     =  cs * a0 + sn * a1;
        T(i, p + 1) = -sn * a0 + cs * a1;

        const Scalar u0 = U(i, p);
        const Scalar u1 = U(i, p + 1);
        U(i, p)     =  cs * u0 + sn * u1;
        U(i, p + 1) = -sn * u0 + cs * u1;
    }

    // Enforce exact zero on the new subdiagonal (1x1/1x1 always succeeds).
    T(p + 1, p) = Scalar{0};
    pivot_ratio_out = Scalar{1};
    return true;
}

// --- Generalized Bai-Demmel swap via 4x4 Kronecker Sylvester ---
//
// Handles 1x1/2x2, 2x2/1x1, and 2x2/2x2 in a unified code path. Builds the
// Kronecker operator K = I_{n2} kron A11 - A22^T kron I_{n1} on a fixed 4x4
// Eigen matrix (padding unused entries with zero); the rank-revealing
// ColPivHouseholderQR pivots unused columns to the back so
// R(m-1, m-1) / R(0, 0) is the correct conditioning signal regardless of
// the active subsystem size m = n1 * n2. Rejection applies both
// conditioning-signal and LAPACK-style post-swap residual tests; when
// accepted, the tentative similarity is written back into T and
// accumulated into U.

template <typename Scalar, int N>
auto swap_real_schur_2x2_general(Eigen::Matrix<Scalar, N, N>& T,
                                 Eigen::Matrix<Scalar, N, N>& U,
                                 int p, int n1, int n2,
                                 Scalar& pivot_ratio_out) -> bool
{
    using std::abs;

    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    const int m = n1 * n2;   // 2, 2, or 4
    const int w = n1 + n2;   // 2, 3, or 4

    // --- Assemble the Kronecker system K * vec(X) = vec(A12) on fixed 4x4. ---
    Eigen::Matrix<Scalar, 4, 4> K   = Eigen::Matrix<Scalar, 4, 4>::Zero();
    Eigen::Matrix<Scalar, 4, 1> rhs = Eigen::Matrix<Scalar, 4, 1>::Zero();

    // vec-by-columns index of X(i, j) is row = j * n1 + i.
    for (int j = 0; j < n2; ++j)
    {
        for (int i = 0; i < n1; ++i)
        {
            const int row = j * n1 + i;
            // (I_{n2} kron A11) * vec(X): row += sum_k A11(i, k) * X(k, j)
            for (int k = 0; k < n1; ++k)
                K(row, j * n1 + k) += T(p + i, p + k);
            // -(A22^T kron I_{n1}) * vec(X): row -= sum_l A22(l, j) * X(i, l)
            for (int l = 0; l < n2; ++l)
                K(row, l * n1 + i) -= T(p + n1 + l, p + n1 + j);
            rhs(row) = T(p + i, p + n1 + j);
        }
    }

    // --- Rank-revealing QR gives the conditioning signal for free. ---
    Eigen::ColPivHouseholderQR<Eigen::Matrix<Scalar, 4, 4>> qr(K);
    const auto R_mat = qr.matrixR();
    const Scalar r_first = abs(R_mat(0,     0));
    const Scalar r_last  = abs(R_mat(m - 1, m - 1));
    pivot_ratio_out = (r_first > Scalar{0}) ? (r_last / r_first) : Scalar{0};

    const Scalar scale_T = T.cwiseAbs().maxCoeff();
    // Structural rejection: Sylvester operator effectively singular.
    if (pivot_ratio_out < eps)
        return false;

    const Eigen::Matrix<Scalar, 4, 1> x_full = qr.solve(rhs);

    // Unpack X (n1 x n2) from the column-stacked solution.
    Eigen::Matrix<Scalar, 2, 2> X = Eigen::Matrix<Scalar, 2, 2>::Zero();
    for (int j = 0; j < n2; ++j)
        for (int i = 0; i < n1; ++i)
            X(i, j) = x_full(j * n1 + i);

    // --- Build orthogonal Q from QR of G = [[-X]; I_{n2}]  (shape (n1+n2) x n2). ---
    Eigen::Matrix<Scalar, 4, 2> G = Eigen::Matrix<Scalar, 4, 2>::Zero();
    for (int j = 0; j < n2; ++j)
    {
        for (int i = 0; i < n1; ++i)
            G(i, j) = -X(i, j);
        G(n1 + j, j) = Scalar{1};
    }

    // Fixed 4x2 Householder QR; reconstruct the (w x w) Q from the full
    // 4x4 householderQ() by truncating to the active rows/cols and padding
    // the trailing block with identity for clean 4x4 fused arithmetic below.
    Eigen::HouseholderQR<Eigen::Matrix<Scalar, 4, 2>> qr_G(G);
    const Eigen::Matrix<Scalar, 4, 4> Q_full = qr_G.householderQ();
    Eigen::Matrix<Scalar, 4, 4> Q = Eigen::Matrix<Scalar, 4, 4>::Identity();
    for (int i = 0; i < w; ++i)
        for (int j = 0; j < w; ++j)
            Q(i, j) = Q_full(i, j);

    // --- Tentative similarity on the w x w block (computed on fixed 4x4). ---
    Eigen::Matrix<Scalar, 4, 4> A_block = Eigen::Matrix<Scalar, 4, 4>::Zero();
    for (int i = 0; i < w; ++i)
        for (int j = 0; j < w; ++j)
            A_block(i, j) = T(p + i, p + j);

    // Compute A_new = Q^T * A_block * Q on the active w x w sub-region
    // without Eigen slice expressions (avoids dynamic-size assignment into
    // fixed 4x4 targets).
    Eigen::Matrix<Scalar, 4, 4> QtA = Eigen::Matrix<Scalar, 4, 4>::Zero();
    for (int i = 0; i < w; ++i)
        for (int j = 0; j < w; ++j)
        {
            Scalar acc = Scalar{0};
            for (int k = 0; k < w; ++k)
                acc += Q(k, i) * A_block(k, j);  // Q^T(i, k) = Q(k, i)
            QtA(i, j) = acc;
        }
    Eigen::Matrix<Scalar, 4, 4> A_new = Eigen::Matrix<Scalar, 4, 4>::Zero();
    for (int i = 0; i < w; ++i)
        for (int j = 0; j < w; ++j)
        {
            Scalar acc = Scalar{0};
            for (int k = 0; k < w; ++k)
                acc += QtA(i, k) * Q(k, j);
            A_new(i, j) = acc;
        }

    // --- Post-swap residual rejection (LAPACK SLAEXC, Bai-Demmel 1993 p.79). ---
    Scalar residual_norm = Scalar{0};
    for (int i = n2; i < w; ++i)
        for (int j = 0; j < n2; ++j)
            residual_norm = std::max(residual_norm, abs(A_new(i, j)));
    const Scalar reject_threshold =
        detail_constants::rejection_multiplier<Scalar>() * eps * scale_T;
    if (residual_norm > reject_threshold)
        return false;

    // --- Accept the swap: write A_new back into T with exact zero on the
    //     new (2,1) block; propagate similarity to off-block rows/cols; and
    //     accumulate the similarity into U. ---

    // Write back the w x w block; zero the new (2,1) block exactly.
    for (int i = 0; i < w; ++i)
        for (int j = 0; j < w; ++j)
            T(p + i, p + j) = A_new(i, j);
    for (int i = n2; i < w; ++i)
        for (int j = 0; j < n2; ++j)
            T(p + i, p + j) = Scalar{0};

    // Rows above/below the w-block (i.e. other cols within rows p..p+w-1
    // were already covered by the block assignment above; the similarity
    // on the full T also needs to update cols outside [p, p+w) for rows
    // [p, p+w) and rows outside [p, p+w) for cols [p, p+w)).

    // Left-multiply rows (p .. p+w-1) by Q^T across all N columns,
    // skipping cols in [p, p+w) which were already handled by the block
    // similarity above.
    for (int j = 0; j < N; ++j)
    {
        if (j >= p && j < p + w)
            continue;  // block region already handled
        Scalar tmp[4] = {Scalar{0}, Scalar{0}, Scalar{0}, Scalar{0}};
        for (int i = 0; i < w; ++i)
            for (int k = 0; k < w; ++k)
                tmp[i] += Q(k, i) * T(p + k, j);  // (Q^T)(i, k) = Q(k, i)
        for (int i = 0; i < w; ++i)
            T(p + i, j) = tmp[i];
    }

    // Right-multiply cols (p .. p+w-1) by Q across all N rows, skipping
    // rows in [p, p+w) already handled by the block similarity.
    for (int i = 0; i < N; ++i)
    {
        if (i >= p && i < p + w)
            continue;
        Scalar tmp[4] = {Scalar{0}, Scalar{0}, Scalar{0}, Scalar{0}};
        for (int j = 0; j < w; ++j)
            for (int k = 0; k < w; ++k)
                tmp[j] += T(i, p + k) * Q(k, j);
        for (int j = 0; j < w; ++j)
            T(i, p + j) = tmp[j];
    }

    // Accumulate Q into U (U_new = U * Q) on cols (p .. p+w-1) across all rows.
    for (int i = 0; i < N; ++i)
    {
        Scalar tmp[4] = {Scalar{0}, Scalar{0}, Scalar{0}, Scalar{0}};
        for (int j = 0; j < w; ++j)
            for (int k = 0; k < w; ++k)
                tmp[j] += U(i, p + k) * Q(k, j);
        for (int j = 0; j < w; ++j)
            U(i, p + j) = tmp[j];
    }

    return true;
}

// --- Swap dispatcher: runtime branch on (n1, n2). ---
//
// The Cond tag is accepted for future policy dispatch. Only
// pivot_ratio_conditioning is reachable today; hager_higham_conditioning is
// blocked by the static_assert inside reorder_real_schur.

template <typename Scalar, int N, typename Cond>
auto swap_real_schur_blocks(Eigen::Matrix<Scalar, N, N>& T,
                            Eigen::Matrix<Scalar, N, N>& U,
                            int p, int n1, int n2,
                            Cond /*tag*/,
                            Scalar& pivot_ratio_out) -> bool
{
    if (n1 == 1 && n2 == 1)
        return swap_real_schur_1x1<Scalar, N>(T, U, p, pivot_ratio_out);
    return swap_real_schur_2x2_general<Scalar, N>(T, U, p, n1, n2, pivot_ratio_out);
}

}

#endif
