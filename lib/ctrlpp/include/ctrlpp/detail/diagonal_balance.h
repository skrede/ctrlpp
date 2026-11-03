#ifndef HPP_GUARD_CTRLPP_DETAIL_DIAGONAL_BALANCE_H
#define HPP_GUARD_CTRLPP_DETAIL_DIAGONAL_BALANCE_H

// LAPACK DGEBAL phase 2 without the permutation phase: equilibrates the off-diagonal row and column
// absolute sums by the diagonal similarity M' = D^-1 M D. Power-of-two factors make the similarity exact, so
// M' has M's eigenvalues.
// @cite lapack_dgebal : Reference LAPACK SRC/dgebal.f phase 2 algorithm

#include <Eigen/Core>

#include <cmath>
#include <limits>
#include <cstdint>
#include <algorithm>

namespace ctrlpp::detail
{

// DGEBAL's SCLFAC, FACTOR, SFMIN2 and SFMAX2.
template <typename Scalar>
struct dgebal_limits
{
    Scalar sclfac;
    Scalar factor;
    Scalar sfmin2;
    Scalar sfmax2;
};

template <typename Scalar>
constexpr dgebal_limits<Scalar> make_dgebal_limits()
{
    const Scalar sclfac = Scalar{2};
    const Scalar factor_num = Scalar{19};
    const Scalar factor_den = Scalar{20};
    const Scalar factor = factor_num / factor_den;
    const Scalar sfmin = std::numeric_limits<Scalar>::min();
    const Scalar sfmax = Scalar{1} / sfmin;
    const Scalar sfmin2 = sfmin * sclfac;
    const Scalar sfmax2 = sfmax / sclfac;
    return {sclfac, factor, sfmin2, sfmax2};
}

// The power of two that brings the column sum c and the row sum r closest, applied to both.
template <typename Scalar>
Scalar dgebal_scale_factor(Scalar& c, Scalar& r)
{
    constexpr dgebal_limits<Scalar> limits = make_dgebal_limits<Scalar>();
    Scalar g = r / limits.sclfac;
    Scalar f = Scalar{1};
    while(c < g && std::max(f, c) < limits.sfmax2 && std::min(r, g) > limits.sfmin2)
    {
        f *= limits.sclfac;
        c *= limits.sclfac;
        r /= limits.sclfac;
        g /= limits.sclfac;
    }
    g = c / limits.sclfac;
    while(g >= r && std::max(r, g) < limits.sfmax2 && std::min({f, c, g}) > limits.sfmin2)
    {
        f /= limits.sclfac;
        c /= limits.sclfac;
        g /= limits.sclfac;
        r *= limits.sclfac;
    }
    return f;
}

template <typename Scalar>
struct off_diagonal_sums
{
    Scalar r;
    Scalar c;
};

template <typename Scalar, int N>
off_diagonal_sums<Scalar> off_diagonal_sums_of(const Eigen::Matrix<Scalar, N, N>& M, Eigen::Index i)
{
    off_diagonal_sums<Scalar> sums{Scalar{0}, Scalar{0}};
    for(Eigen::Index j = 0; j < N; ++j)
    {
        if(j == i)
            continue;
        sums.r += std::abs(M(i, j));
        sums.c += std::abs(M(j, i));
    }
    return sums;
}

// The sweep cap is a safety ceiling, not a prediction: the count of power-of-two scalings an index can take
// inside the representable range, once per index. Stopping at it is safe because every applied step updates M
// and D_out together, so M == D^-1 M_original D holds exactly at any cut point.
template <typename Scalar, int N>
void balance_diagonal(Eigen::Matrix<Scalar, N, N>& M, Eigen::Matrix<Scalar, N, 1>& D_out)
{
    constexpr std::int32_t max_sweeps = N * (std::numeric_limits<Scalar>::max_exponent - std::numeric_limits<Scalar>::min_exponent);
    constexpr dgebal_limits<Scalar> limits = make_dgebal_limits<Scalar>();

    D_out = Eigen::Matrix<Scalar, N, 1>::Ones();
    bool noconv = true;
    for(std::int32_t sweep = 0; noconv && sweep < max_sweeps; ++sweep)
    {
        noconv = false;
        for(Eigen::Index i = 0; i < N; ++i)
        {
            auto [r, c] = off_diagonal_sums_of<Scalar, N>(M, i);
            if(r == Scalar{0} || c == Scalar{0})
                continue;
            const Scalar s = c + r;
            const Scalar f = dgebal_scale_factor(c, r);
            if((c + r) < limits.factor * s)
            {
                noconv = true;
                D_out(i) *= f;
                M.row(i) /= f;
                M.col(i) *= f;
            }
        }
    }
}
}

#endif
