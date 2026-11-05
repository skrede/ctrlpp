#ifndef HPP_GUARD_CTRLPP_DETAIL_SCHUR_REORDER_H
#define HPP_GUARD_CTRLPP_DETAIL_SCHUR_REORDER_H

// Real-Schur block reordering by the Bai and Demmel predicated swap of adjacent diagonal blocks (LAPACK
// DTREXC, DLAEXC), with DLANV2 Murnaghan standardization of every 2x2 block a swap touches.
// @cite bai_demmel_1993 -- Bai & Demmel, "On swapping diagonal blocks in real Schur form", 1993
// @cite laub1979 -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979

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

// Rejects an ill-conditioned swap by |R(last, last)| / |R(0, 0)| of the rank-revealing QR of the Sylvester
// operator, which the swap computes anyway.
struct pivot_ratio_conditioning
{
};

// Hager-Higham 1-norm estimate of the inverse Sylvester operator; not implemented, so selecting it fails to
// compile.
struct hager_higham_conditioning
{
};

// The solvers constrain their trailing policy parameter on this, so overload resolution rejects a foreign
// type in that slot instead of finding an ambiguity with the overloads taking one more matrix.
template <typename T>
concept conditioning_policy = std::same_as<T, pivot_ratio_conditioning> || std::same_as<T, hager_higham_conditioning>;

namespace detail_constants
{

// LAPACK DLANV2 MULTPL: the real-versus-complex-pair decision multiplier.
template <typename Scalar>
constexpr Scalar dlanv2_multpl() noexcept
{
    return Scalar{4};
}

// LAPACK DLAEXC TEN: the post-swap residual threshold multiplier.
template <typename Scalar>
constexpr Scalar rejection_multiplier() noexcept
{
    return Scalar{10};
}

}

// `placed` counts the selected eigenvalues moved to the top; `complete` is false once a swap was rejected;
// `subspace_separation` is the smallest pivot ratio over the swaps, the analogue of LAPACK's SEP.
template <typename Scalar>
struct reorder_result
{
    reorder_result()
        : complete{true}
        , placed{0}
        , subspace_separation{Scalar{1}}
    {
    }

    bool complete;
    Eigen::Index placed;
    Scalar subspace_separation;
};

template <typename Scalar, int N>
void standardize_2x2_block(Eigen::Matrix<Scalar, N, N>& T, Eigen::Matrix<Scalar, N, N>& U, Eigen::Index p);

template <typename Scalar, int N, typename Cond>
bool swap_real_schur_blocks(Eigen::Matrix<Scalar, N, N>& T, Eigen::Matrix<Scalar, N, N>& U, Eigen::Index p, Eigen::Index n1, Eigen::Index n2, Cond, Scalar& pivot_ratio_out);

template <typename Scalar, int N, typename Predicate, typename Cond = pivot_ratio_conditioning>
reorder_result<Scalar> reorder_real_schur(Eigen::Matrix<Scalar, N, N>& T, Eigen::Matrix<Scalar, N, N>& U, Predicate&& predicate, Cond tag = {})
{
    static_assert(
            std::is_same_v<Cond, pivot_ratio_conditioning> || std::is_same_v<Cond, hager_higham_conditioning>, "Cond must be pivot_ratio_conditioning or hager_higham_conditioning");
    if constexpr(std::is_same_v<Cond, hager_higham_conditioning>)
    {
        static_assert(!std::is_same_v<Cond, hager_higham_conditioning>, "hager_higham_conditioning is not implemented; use pivot_ratio_conditioning");
    }

    reorder_result<Scalar> r;

    const Scalar scale_T = T.cwiseAbs().maxCoeff();

    auto block_size_at = [&](Eigen::Index pos) -> Eigen::Index { return quasi_triangular_block_size<Scalar, N>(T, pos, scale_T); };

    // `first` is the larger real root of a real pair and shares its modulus and real part with the other
    // root of a conjugate pair, so the predicate asks about `first` alone.
    auto block_matches = [&](Eigen::Index pos, Eigen::Index n_block) -> bool { return predicate(quasi_triangular_block_spectrum<Scalar, N>(T, pos, n_block).first); };

    Eigen::Index placed_pos = 0;
    while(placed_pos < N)
    {
        Eigen::Index scan = placed_pos;
        Eigen::Index scan_nb = 1;
        bool found = false;
        while(scan < N)
        {
            const Eigen::Index nb = block_size_at(scan);
            if(block_matches(scan, nb))
            {
                scan_nb = nb;
                found = true;
                break;
            }
            scan += nb;
        }

        if(!found)
            break;

        Eigen::Index cur = scan;
        while(cur > placed_pos)
        {
            // A 2x2 block ending at cur - 1 starts at cur - 2.
            const Eigen::Index left_nb = (cur >= 2) ? block_size_at(cur - 2) : 1;
            const Eigen::Index left_pos = cur - left_nb;

            Scalar pivot_ratio = Scalar{1};
            const bool ok = swap_real_schur_blocks<Scalar, N, Cond>(T, U, left_pos, left_nb, scan_nb, tag, pivot_ratio);

            if(!ok)
            {
                r.complete = false;
                r.subspace_separation = std::min(r.subspace_separation, pivot_ratio);
                break;
            }

            r.subspace_separation = std::min(r.subspace_separation, pivot_ratio);

            if(scan_nb == 2)
                standardize_2x2_block<Scalar, N>(T, U, left_pos);
            if(left_nb == 2)
                standardize_2x2_block<Scalar, N>(T, U, left_pos + scan_nb);

            cur = left_pos;
            // Standardization can split a 2x2 block whose conjugate pair degenerated into two real
            // eigenvalues, and the split block may no longer match.
            scan_nb = block_size_at(cur);
            if(!block_matches(cur, scan_nb))
                break;
        }

        if(cur == placed_pos)
        {
            r.placed += scan_nb;
            placed_pos += scan_nb;
        }
        else
        {
            // Resuming past the stuck block guarantees termination.
            placed_pos = cur + block_size_at(cur);
        }
    }

    return r;
}

// LAPACK DLANV2: a real pair becomes upper triangular; a conjugate pair gets equal diagonals and
// opposite-sign off-diagonals. The rotation is accumulated into U.
template <typename Scalar, int N>
void standardize_2x2_block(Eigen::Matrix<Scalar, N, N>& T, Eigen::Matrix<Scalar, N, N>& U, Eigen::Index p)
{
    using std::abs;
    using std::sqrt;
    using std::hypot;
    using std::copysign;

    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    const Scalar multpl = detail_constants::dlanv2_multpl<Scalar>();

    const Scalar a = T(p, p);
    const Scalar b = T(p, p + 1);
    const Scalar c = T(p + 1, p);
    const Scalar d = T(p + 1, p + 1);

    Scalar cs = Scalar{1};
    Scalar sn = Scalar{0};
    bool enforce_equal_diag = false;
    Scalar ad_mean = Scalar{0};

    if(c == Scalar{0})
    {
        cs = Scalar{1};
        sn = Scalar{0};
    }
    else if(b == Scalar{0})
    {
        cs = Scalar{0};
        sn = Scalar{1};
    }
    else if(a == d && ((b > Scalar{0}) != (c > Scalar{0})))
    {
        cs = Scalar{1};
        sn = Scalar{0};
    }
    else
    {
        const Scalar temp = a - d;
        const Scalar p_val = temp / Scalar{2};
        const Scalar bcmax = std::max(abs(b), abs(c));
        const Scalar sgn_b = (b < Scalar{0}) ? Scalar{-1} : Scalar{1};
        const Scalar sgn_c = (c < Scalar{0}) ? Scalar{-1} : Scalar{1};
        const Scalar bcmis = std::min(abs(b), abs(c)) * sgn_b * sgn_c;
        const Scalar scl = std::max(abs(p_val), bcmax);
        const Scalar z_val = (p_val / scl) * p_val + (bcmax / scl) * bcmis;

        if(z_val >= multpl * eps)
        {
            Scalar z = p_val + copysign(sqrt(scl) * sqrt(z_val), p_val);
            Scalar tau = hypot(c, z);
            cs = z / tau;
            sn = c / tau;
        }
        else
        {
            const Scalar sigma = b + c;
            const Scalar sgn_sig = (sigma < Scalar{0}) ? Scalar{-1} : Scalar{1};
            Scalar tau = hypot(sigma, temp);
            cs = sqrt((Scalar{1} + abs(sigma) / tau) / Scalar{2});
            sn = -(p_val / (tau * cs)) * sgn_sig;
            enforce_equal_diag = true;
            ad_mean = (a + d) / Scalar{2};
        }
    }

    for(Eigen::Index j = 0; j < N; ++j)
    {
        const Scalar t0 = T(p, j);
        const Scalar t1 = T(p + 1, j);
        T(p, j) = cs * t0 + sn * t1;
        T(p + 1, j) = -sn * t0 + cs * t1;
    }
    for(Eigen::Index i = 0; i < N; ++i)
    {
        const Scalar t0 = T(i, p);
        const Scalar t1 = T(i, p + 1);
        T(i, p) = cs * t0 + sn * t1;
        T(i, p + 1) = -sn * t0 + cs * t1;

        const Scalar u0 = U(i, p);
        const Scalar u1 = U(i, p + 1);
        U(i, p) = cs * u0 + sn * u1;
        U(i, p + 1) = -sn * u0 + cs * u1;
    }

    // Exactly equal diagonals suppress the drift of the sigma-derived rotation.
    if(enforce_equal_diag)
    {
        T(p, p) = ad_mean;
        T(p + 1, p + 1) = ad_mean;
    }
}

// The 1x1/1x1 swap always succeeds (Bai and Demmel, p. 79). Coinciding eigenvalues make it a no-op, as in
// LAPACK; otherwise (cs, sn) = (t12, -diff) / hypot(t12, diff) zeroes the new subdiagonal without overflow
// (Golub and Van Loan, Alg. 5.1.3).
template <typename Scalar, int N>
bool swap_real_schur_1x1(Eigen::Matrix<Scalar, N, N>& T, Eigen::Matrix<Scalar, N, N>& U, Eigen::Index p, Scalar& pivot_ratio_out)
{
    using std::abs;
    using std::sqrt;
    using std::hypot;

    const Scalar eps = std::numeric_limits<Scalar>::epsilon();

    const Scalar t11 = T(p, p);
    const Scalar t22 = T(p + 1, p + 1);
    const Scalar t12 = T(p, p + 1);
    const Scalar diff = t11 - t22;

    if(abs(diff) < eps * std::max(abs(t11), abs(t22)))
    {
        pivot_ratio_out = Scalar{1};
        return true;
    }

    const Scalar r = hypot(t12, diff);
    const Scalar cs = t12 / r;
    const Scalar sn = -diff / r;

    for(Eigen::Index j = 0; j < N; ++j)
    {
        const Scalar a0 = T(p, j);
        const Scalar a1 = T(p + 1, j);
        T(p, j) = cs * a0 + sn * a1;
        T(p + 1, j) = -sn * a0 + cs * a1;
    }
    for(Eigen::Index i = 0; i < N; ++i)
    {
        const Scalar a0 = T(i, p);
        const Scalar a1 = T(i, p + 1);
        T(i, p) = cs * a0 + sn * a1;
        T(i, p + 1) = -sn * a0 + cs * a1;

        const Scalar u0 = U(i, p);
        const Scalar u1 = U(i, p + 1);
        U(i, p) = cs * u0 + sn * u1;
        U(i, p + 1) = -sn * u0 + cs * u1;
    }

    T(p + 1, p) = Scalar{0};
    pivot_ratio_out = Scalar{1};
    return true;
}

// The 1x2, 2x1 and 2x2 swaps share one path: the Sylvester operator K = I kron A11 - A22^T kron I sits in a
// zero-padded fixed 4x4, and column pivoting moves the unused columns to the back, so R(m - 1, m - 1) / R(0, 0)
// is the conditioning signal for every active size m = n1 * n2. The swap is rejected on that signal or on the
// post-swap residual (LAPACK DLAEXC, Bai and Demmel, p. 79).
template <typename Scalar, int N>
bool swap_real_schur_2x2_general(Eigen::Matrix<Scalar, N, N>& T, Eigen::Matrix<Scalar, N, N>& U, Eigen::Index p, Eigen::Index n1, Eigen::Index n2, Scalar& pivot_ratio_out)
{
    using std::abs;

    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    const Eigen::Index m = n1 * n2;
    const Eigen::Index w = n1 + n2;

    Eigen::Matrix<Scalar, 4, 4> K = Eigen::Matrix<Scalar, 4, 4>::Zero();
    Eigen::Matrix<Scalar, 4, 1> rhs = Eigen::Matrix<Scalar, 4, 1>::Zero();

    // vec(X) stacks columns: X(i, j) is entry j * n1 + i.
    for(Eigen::Index j = 0; j < n2; ++j)
    {
        for(Eigen::Index i = 0; i < n1; ++i)
        {
            const Eigen::Index row = j * n1 + i;
            for(Eigen::Index k = 0; k < n1; ++k)
                K(row, j * n1 + k) += T(p + i, p + k);
            for(Eigen::Index l = 0; l < n2; ++l)
                K(row, l * n1 + i) -= T(p + n1 + l, p + n1 + j);
            rhs(row) = T(p + i, p + n1 + j);
        }
    }

    Eigen::ColPivHouseholderQR<Eigen::Matrix<Scalar, 4, 4>> qr(K);
    const auto R_mat = qr.matrixR();
    const Scalar r_first = abs(R_mat(0, 0));
    const Scalar r_last = abs(R_mat(m - 1, m - 1));
    pivot_ratio_out = (r_first > Scalar{0}) ? (r_last / r_first) : Scalar{0};

    const Scalar scale_T = T.cwiseAbs().maxCoeff();
    if(pivot_ratio_out < eps)
        return false;

    const Eigen::Matrix<Scalar, 4, 1> x_full = qr.solve(rhs);

    Eigen::Matrix<Scalar, 2, 2> X = Eigen::Matrix<Scalar, 2, 2>::Zero();
    for(Eigen::Index j = 0; j < n2; ++j)
        for(Eigen::Index i = 0; i < n1; ++i)
            X(i, j) = x_full(j * n1 + i);

    Eigen::Matrix<Scalar, 4, 2> G = Eigen::Matrix<Scalar, 4, 2>::Zero();
    for(Eigen::Index j = 0; j < n2; ++j)
    {
        for(Eigen::Index i = 0; i < n1; ++i)
            G(i, j) = -X(i, j);
        G(n1 + j, j) = Scalar{1};
    }

    // Q is the active w x w corner of the Householder Q of [-X; I], padded with the identity.
    Eigen::HouseholderQR<Eigen::Matrix<Scalar, 4, 2>> qr_G(G);
    const Eigen::Matrix<Scalar, 4, 4> Q_full = qr_G.householderQ();
    Eigen::Matrix<Scalar, 4, 4> Q = Eigen::Matrix<Scalar, 4, 4>::Identity();
    for(Eigen::Index i = 0; i < w; ++i)
        for(Eigen::Index j = 0; j < w; ++j)
            Q(i, j) = Q_full(i, j);

    Eigen::Matrix<Scalar, 4, 4> A_block = Eigen::Matrix<Scalar, 4, 4>::Zero();
    for(Eigen::Index i = 0; i < w; ++i)
        for(Eigen::Index j = 0; j < w; ++j)
            A_block(i, j) = T(p + i, p + j);

    // Explicit loops rather than slice expressions avoid a dynamic-size assignment into a fixed 4x4.
    Eigen::Matrix<Scalar, 4, 4> QtA = Eigen::Matrix<Scalar, 4, 4>::Zero();
    for(Eigen::Index i = 0; i < w; ++i)
        for(Eigen::Index j = 0; j < w; ++j)
        {
            Scalar acc = Scalar{0};
            for(Eigen::Index k = 0; k < w; ++k)
                acc += Q(k, i) * A_block(k, j);
            QtA(i, j) = acc;
        }
    Eigen::Matrix<Scalar, 4, 4> A_new = Eigen::Matrix<Scalar, 4, 4>::Zero();
    for(Eigen::Index i = 0; i < w; ++i)
        for(Eigen::Index j = 0; j < w; ++j)
        {
            Scalar acc = Scalar{0};
            for(Eigen::Index k = 0; k < w; ++k)
                acc += QtA(i, k) * Q(k, j);
            A_new(i, j) = acc;
        }

    Scalar residual_norm = Scalar{0};
    for(Eigen::Index i = n2; i < w; ++i)
        for(Eigen::Index j = 0; j < n2; ++j)
            residual_norm = std::max(residual_norm, abs(A_new(i, j)));
    const Scalar reject_threshold = detail_constants::rejection_multiplier<Scalar>() * eps * scale_T;
    if(residual_norm > reject_threshold)
        return false;

    for(Eigen::Index i = 0; i < w; ++i)
        for(Eigen::Index j = 0; j < w; ++j)
            T(p + i, p + j) = A_new(i, j);
    for(Eigen::Index i = n2; i < w; ++i)
        for(Eigen::Index j = 0; j < n2; ++j)
            T(p + i, p + j) = Scalar{0};

    for(Eigen::Index j = 0; j < N; ++j)
    {
        if(j >= p && j < p + w)
            continue;
        Eigen::Matrix<Scalar, 4, 1> tmp = Eigen::Matrix<Scalar, 4, 1>::Zero();
        for(Eigen::Index i = 0; i < w; ++i)
            for(Eigen::Index k = 0; k < w; ++k)
                tmp(i) += Q(k, i) * T(p + k, j);
        for(Eigen::Index i = 0; i < w; ++i)
            T(p + i, j) = tmp(i);
    }

    for(Eigen::Index i = 0; i < N; ++i)
    {
        if(i >= p && i < p + w)
            continue;
        Eigen::Matrix<Scalar, 4, 1> tmp = Eigen::Matrix<Scalar, 4, 1>::Zero();
        for(Eigen::Index j = 0; j < w; ++j)
            for(Eigen::Index k = 0; k < w; ++k)
                tmp(j) += T(i, p + k) * Q(k, j);
        for(Eigen::Index j = 0; j < w; ++j)
            T(i, p + j) = tmp(j);
    }

    for(Eigen::Index i = 0; i < N; ++i)
    {
        Eigen::Matrix<Scalar, 4, 1> tmp = Eigen::Matrix<Scalar, 4, 1>::Zero();
        for(Eigen::Index j = 0; j < w; ++j)
            for(Eigen::Index k = 0; k < w; ++k)
                tmp(j) += U(i, p + k) * Q(k, j);
        for(Eigen::Index j = 0; j < w; ++j)
            U(i, p + j) = tmp(j);
    }

    return true;
}

template <typename Scalar, int N, typename Cond>
bool swap_real_schur_blocks(Eigen::Matrix<Scalar, N, N>& T, Eigen::Matrix<Scalar, N, N>& U, Eigen::Index p, Eigen::Index n1, Eigen::Index n2, Cond, Scalar& pivot_ratio_out)
{
    if(n1 == 1 && n2 == 1)
        return swap_real_schur_1x1<Scalar, N>(T, U, p, pivot_ratio_out);
    return swap_real_schur_2x2_general<Scalar, N>(T, U, p, n1, n2, pivot_ratio_out);
}

}

#endif
