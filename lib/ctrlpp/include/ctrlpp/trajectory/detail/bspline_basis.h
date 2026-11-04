#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_BSPLINE_BASIS_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_BSPLINE_BASIS_H

#include <Eigen/Dense>

#include <array>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar>
inline std::int32_t find_span(Scalar t, std::vector<Scalar> const& U, std::int32_t n, std::int32_t degree)
{
    if(t >= U[static_cast<std::size_t>(n + 1)])
        return n;

    auto it = std::upper_bound(U.begin() + degree, U.begin() + n + 2, t);
    auto k = static_cast<std::int32_t>(it - U.begin()) - 1;
    k = std::clamp(k, degree, n);
    return k;
}

// @cite biagiotti2009 -- Sec. 4.5
template <int Degree, typename Scalar>
inline Scalar de_boor(Scalar t, std::vector<Scalar> const& points, std::vector<Scalar> const& U, std::int32_t degree)
{
    auto const n = static_cast<std::int32_t>(points.size()) - 1;
    std::int32_t k = find_span(t, U, n, degree);

    // Sized by the class degree, which bounds every call's degree, so evaluation does not allocate.
    std::array<Scalar, static_cast<std::size_t>(Degree) + 1> d{};
    for(std::int32_t j = 0; j <= degree; ++j)
        d[static_cast<std::size_t>(j)] = points[static_cast<std::size_t>(k - degree + j)];

    for(std::int32_t r = 1; r <= degree; ++r)
    {
        for(std::int32_t j = degree; j >= r; --j)
        {
            auto const idx = k - degree + j;
            auto const denom = U[static_cast<std::size_t>(idx + degree - r + 1)] - U[static_cast<std::size_t>(idx)];
            // 0/0 := 0, Piegl & Tiller Eq. 2.5; only t == U[n + 1] with U[n] == U[n + 1] reaches it.
            auto const alpha = denom == Scalar{0} ? Scalar{0} : (t - U[static_cast<std::size_t>(idx)]) / denom;
            d[static_cast<std::size_t>(j)] = (Scalar{1} - alpha) * d[static_cast<std::size_t>(j - 1)] + alpha * d[static_cast<std::size_t>(j)];
        }
    }

    return d[static_cast<std::size_t>(degree)];
}

template <typename Scalar>
Scalar bspline_basis_degree_zero(std::int32_t i, Scalar t, std::vector<Scalar> const& U)
{
    auto const ui = U[static_cast<std::size_t>(i)];
    auto const ui1 = U[static_cast<std::size_t>(i + 1)];
    // At the right end of the domain a span of non-zero width also takes its right boundary (clamped endpoint rule).
    auto const t_max = U.back();
    if(t >= ui && t < ui1)
        return Scalar{1};
    if(t == t_max && ui < ui1 && t >= ui && t <= ui1)
        return Scalar{1};
    return Scalar{0};
}

// Cox-de Boor recursion.
// @cite biagiotti2009 -- Sec. 4.5
template <typename Scalar>
Scalar bspline_basis(std::int32_t i, std::int32_t p, Scalar t, std::vector<Scalar> const& U)
{
    if(p == 0)
        return bspline_basis_degree_zero(i, t, U);

    Scalar left{0};
    Scalar right{0};

    auto const denom_left = U[static_cast<std::size_t>(i + p)] - U[static_cast<std::size_t>(i)];
    // 0/0 := 0, Piegl & Tiller Eq. 2.5; outside the lower basis' support the ratio may overflow, so skip it.
    if(auto const lower = bspline_basis(i, p - 1, t, U); denom_left != Scalar{0} && lower != Scalar{0})
        left = (t - U[static_cast<std::size_t>(i)]) / denom_left * lower;

    auto const denom_right = U[static_cast<std::size_t>(i + p + 1)] - U[static_cast<std::size_t>(i + 1)];
    if(auto const lower = bspline_basis(i + 1, p - 1, t, U); denom_right != Scalar{0} && lower != Scalar{0})
        right = (U[static_cast<std::size_t>(i + p + 1)] - t) / denom_right * lower;

    return left + right;
}

// Clamped knots: p + 1 copies of each end time, interior knots by de Boor's averaging
// u_{j+p} = (1/p) sum_{i=j}^{j+p-1} t_i.
template <int Degree, typename Scalar>
inline std::vector<Scalar> averaged_knots(std::vector<Scalar> const& times)
{
    auto const n = static_cast<std::int32_t>(times.size()) - 1;
    constexpr std::int32_t p = Degree;
    auto const m = n + p + 1;
    std::vector<Scalar> knots(static_cast<std::size_t>(m + 1));

    for(std::int32_t i = 0; i <= p; ++i)
        knots[static_cast<std::size_t>(i)] = times[0];
    for(std::int32_t i = m - p; i <= m; ++i)
        knots[static_cast<std::size_t>(i)] = times[static_cast<std::size_t>(n)];

    for(std::int32_t j = 1; j <= n - p; ++j)
    {
        Scalar sum{0};
        for(std::int32_t i = j; i < j + p; ++i)
            sum += times[static_cast<std::size_t>(i)];
        knots[static_cast<std::size_t>(j + p)] = sum / static_cast<Scalar>(p);
    }
    return knots;
}

// Collocation matrix N(i, j) = B_{j,p}(t_i) of the system N P = Q.
template <int Degree, typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> collocation_matrix(std::vector<Scalar> const& times, std::vector<Scalar> const& knots)
{
    auto const n = static_cast<std::int32_t>(times.size()) - 1;
    constexpr std::int32_t p = Degree;
    Eigen::Index const num_basis = n + 1;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> N(num_basis, num_basis);
    N.setZero();

    for(Eigen::Index i = 0; i <= n; ++i)
    {
        auto const t = times[static_cast<std::size_t>(i)];
        for(Eigen::Index j = 0; j <= n; ++j)
            N(i, j) = bspline_basis(static_cast<std::int32_t>(j), p, t, knots);
    }
    return N;
}

}

#endif
