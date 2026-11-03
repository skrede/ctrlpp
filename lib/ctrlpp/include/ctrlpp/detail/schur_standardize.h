#ifndef HPP_GUARD_CTRLPP_DETAIL_SCHUR_STANDARDIZE_H
#define HPP_GUARD_CTRLPP_DETAIL_SCHUR_STANDARDIZE_H

#include "ctrlpp/detail/schur_conditioning.h"

#include <Eigen/Core>

#include <cmath>
#include <limits>
#include <algorithm>

namespace ctrlpp::detail
{

// `enforce_equal_diag` is set on a conjugate pair, whose diagonals are then written as `ad_mean` exactly to
// suppress the drift of the sigma-derived rotation.
template <typename Scalar>
struct dlanv2_rotation
{
    Scalar cs;
    Scalar sn;
    bool enforce_equal_diag;
    Scalar ad_mean;
};

// The rotation G = [[cs, -sn], [sn, cs]] applied as the similarity G^T T G on rows and columns (p, p + 1),
// and accumulated into U.
template <typename Scalar, int N>
void apply_givens_similarity(Eigen::Matrix<Scalar, N, N>& T, Eigen::Matrix<Scalar, N, N>& U, Eigen::Index p, Scalar cs, Scalar sn)
{
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
}

template <typename Scalar>
dlanv2_rotation<Scalar> dlanv2_real_pair_rotation(Scalar c, Scalar p_val, Scalar scl, Scalar z_val)
{
    using std::sqrt;
    using std::hypot;
    using std::copysign;

    Scalar z = p_val + copysign(sqrt(scl) * sqrt(z_val), p_val);
    Scalar tau = hypot(c, z);
    return {z / tau, c / tau, false, Scalar{0}};
}

template <typename Scalar>
dlanv2_rotation<Scalar> dlanv2_conjugate_pair_rotation(Scalar a, Scalar b, Scalar c, Scalar d, Scalar temp, Scalar p_val)
{
    using std::abs;
    using std::sqrt;
    using std::hypot;

    const Scalar sigma = b + c;
    const Scalar sgn_sig = (sigma < Scalar{0}) ? Scalar{-1} : Scalar{1};
    Scalar tau = hypot(sigma, temp);
    const Scalar cs = sqrt((Scalar{1} + abs(sigma) / tau) / Scalar{2});
    const Scalar sn = -(p_val / (tau * cs)) * sgn_sig;
    return {cs, sn, true, (a + d) / Scalar{2}};
}

template <typename Scalar>
dlanv2_rotation<Scalar> dlanv2_general_rotation(Scalar a, Scalar b, Scalar c, Scalar d)
{
    using std::abs;

    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    const Scalar multpl = detail_constants::dlanv2_multpl<Scalar>();

    const Scalar temp = a - d;
    const Scalar p_val = temp / Scalar{2};
    const Scalar bcmax = std::max(abs(b), abs(c));
    const Scalar sgn_b = (b < Scalar{0}) ? Scalar{-1} : Scalar{1};
    const Scalar sgn_c = (c < Scalar{0}) ? Scalar{-1} : Scalar{1};
    const Scalar bcmis = std::min(abs(b), abs(c)) * sgn_b * sgn_c;
    const Scalar scl = std::max(abs(p_val), bcmax);
    const Scalar z_val = (p_val / scl) * p_val + (bcmax / scl) * bcmis;

    if(z_val >= multpl * eps)
        return dlanv2_real_pair_rotation(c, p_val, scl, z_val);
    return dlanv2_conjugate_pair_rotation(a, b, c, d, temp, p_val);
}

// A zero subdiagonal is already triangular, a zero superdiagonal is swapped into it, and equal diagonals
// with off-diagonals of opposite sign are already standardized.
template <typename Scalar>
dlanv2_rotation<Scalar> dlanv2_rotation_of(Scalar a, Scalar b, Scalar c, Scalar d)
{
    if(c == Scalar{0})
        return {Scalar{1}, Scalar{0}, false, Scalar{0}};
    if(b == Scalar{0})
        return {Scalar{0}, Scalar{1}, false, Scalar{0}};
    if(a == d && ((b > Scalar{0}) != (c > Scalar{0})))
        return {Scalar{1}, Scalar{0}, false, Scalar{0}};
    return dlanv2_general_rotation(a, b, c, d);
}

// LAPACK DLANV2: a real pair becomes upper triangular; a conjugate pair gets equal diagonals and
// opposite-sign off-diagonals. The rotation is accumulated into U.
template <typename Scalar, int N>
void standardize_2x2_block(Eigen::Matrix<Scalar, N, N>& T, Eigen::Matrix<Scalar, N, N>& U, Eigen::Index p)
{
    const dlanv2_rotation<Scalar> rotation = dlanv2_rotation_of(T(p, p), T(p, p + 1), T(p + 1, p), T(p + 1, p + 1));
    apply_givens_similarity<Scalar, N>(T, U, p, rotation.cs, rotation.sn);
    if(rotation.enforce_equal_diag)
    {
        T(p, p) = rotation.ad_mean;
        T(p + 1, p + 1) = rotation.ad_mean;
    }
}

}

#endif
