#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_SMOOTHING_SPLINE_SOLVE_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_SMOOTHING_SPLINE_SOLVE_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
// Automatic Machines and Robots", 2009, Sec. 4.4.5

#include "ctrlpp/trajectory/detail/spline_checks.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <limits>
#include <vector>
#include <cstddef>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar>
using smoothing_matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

template <typename Scalar>
using smoothing_vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

template <typename Scalar>
inline Scalar smoothing_weight(Scalar mu)
{
    return Scalar{2} * (Scalar{1} - mu) / (Scalar{3} * mu);
}

// R is diagonally dominant and Q^T Q a Gram matrix, so neither exceeds its largest diagonal entry; Q^T Q's diagonal at
// interior knot j is 1/h_j^2 + (1/h_j + 1/h_{j+1})^2 + 1/h_{j+1}^2. An infinite weight fails the bound like a large one.
template <typename Scalar>
inline bool regularized_system_representable(std::vector<Scalar> const& times, Scalar mu)
{
    auto const n_pts = times.size();
    if(n_pts < 3)
    {
        return true;
    }

    auto const lambda = smoothing_weight(mu);
    auto const bound = std::sqrt(std::numeric_limits<Scalar>::max());

    Scalar largest_diagonal{0};
    for(std::size_t j = 0; j + 2 < n_pts; ++j)
    {
        auto const first = times[j + 1] - times[j];
        auto const second = times[j + 2] - times[j + 1];
        auto const inv_first = Scalar{1} / first;
        auto const inv_second = Scalar{1} / second;
        auto const middle = inv_first + inv_second;

        auto const gram = inv_first * inv_first + middle * middle + inv_second * inv_second;
        auto const smoothness = (first + second) / Scalar{3};
        largest_diagonal = std::max(largest_diagonal, smoothness + lambda * gram);
    }
    return largest_diagonal <= bound;
}

template <typename Scalar>
inline smoothing_matrix<Scalar> smoothness_matrix(std::vector<Scalar> const& h, Eigen::Index mi)
{
    smoothing_matrix<Scalar> R = smoothing_matrix<Scalar>::Zero(mi, mi);
    for(Eigen::Index i = 0; i < mi; ++i)
    {
        auto const ii = static_cast<std::size_t>(i);
        R(i, i) = (h[ii] + h[ii + 1]) / Scalar{3};
        if(i + 1 < mi)
        {
            R(i, i + 1) = h[ii + 1] / Scalar{6};
            R(i + 1, i) = h[ii + 1] / Scalar{6};
        }
    }
    return R;
}

template <typename Scalar>
inline smoothing_matrix<Scalar> second_difference_matrix(std::vector<Scalar> const& h, Eigen::Index ni, Eigen::Index mi)
{
    smoothing_matrix<Scalar> Q = smoothing_matrix<Scalar>::Zero(ni, mi);
    for(Eigen::Index j = 0; j < mi; ++j)
    {
        auto const jj = static_cast<std::size_t>(j);
        Q(j, j) = Scalar{1} / h[jj];
        Q(j + 1, j) = -(Scalar{1} / h[jj] + Scalar{1} / h[jj + 1]);
        Q(j + 2, j) = Scalar{1} / h[jj + 1];
    }
    return Q;
}

template <typename Scalar>
inline smoothing_vector<Scalar> knot_vector(std::vector<Scalar> const& pos)
{
    auto const ni = static_cast<Eigen::Index>(pos.size());
    smoothing_vector<Scalar> q(ni);
    for(Eigen::Index i = 0; i < ni; ++i)
    {
        q(i) = pos[static_cast<std::size_t>(i)];
    }
    return q;
}

// The numerator is a difference of second derivatives that constant curvature cancels to a residual, so a subnormal cubic
// coefficient is refused only when its term reaches past the position's resolution. The reach is ordered so the span is
// never squared alone, which at the largest admitted span overflows.
template <typename Scalar>
inline bool smoothing_coefficients_representable(std::array<Scalar, 4> const& c, Scalar cubic_numerator, Scalar span, Scalar scale)
{
    auto const reach = (std::abs(cubic_numerator) * span / Scalar{6}) * span;
    auto const resolution = spline_position_resolution(scale);

    bool const finite = std::isfinite(c[0]) && std::isfinite(c[1]) && std::isfinite(c[2]) && std::isfinite(c[3]);
    bool const cubic_survived = (reach <= resolution) || std::isnormal(c[3]);
    return finite && cubic_survived;
}

template <typename Scalar>
inline bool build_linear_coefficients(std::vector<std::array<Scalar, 4>>& coeffs, std::vector<Scalar> const& pos, std::vector<Scalar> const& h)
{
    auto const slope = (pos[1] - pos[0]) / h[0];
    coeffs.resize(1);
    coeffs[0][0] = pos[0];
    coeffs[0][1] = slope;
    coeffs[0][2] = Scalar{0};
    coeffs[0][3] = Scalar{0};

    // The curvature terms are assigned rather than divided, so only the slope can leave the type here.
    return std::isfinite(pos[0]) && std::isfinite(slope);
}

// The moment form of the cubic spline coefficients.
// @cite biagiotti2009 -- Sec. 4.4, eq. (4.10)-(4.11)
template <typename Scalar>
inline bool build_smoothing_coefficients(std::vector<std::array<Scalar, 4>>& coeffs, std::vector<Scalar> const& pos, std::vector<Scalar> const& h, std::vector<Scalar> const& d,
        std::vector<Scalar> const& s, std::size_t n)
{
    auto const position_scale = largest_magnitude(pos);
    bool representable = true;
    coeffs.resize(n);
    for(std::size_t i = 0; i < n; ++i)
    {
        auto const si = s[i];
        auto const si1 = s[i + 1];
        auto const di = d[i];
        auto const di1 = d[i + 1];
        auto const hi = h[i];
        auto const cubic_numerator = di1 - di;

        coeffs[i][0] = si;
        coeffs[i][1] = (si1 - si) / hi - hi * (Scalar{2} * di + di1) / Scalar{6};
        coeffs[i][2] = di / Scalar{2};
        coeffs[i][3] = cubic_numerator / (Scalar{6} * hi);

        representable = representable && smoothing_coefficients_representable(coeffs[i], cubic_numerator, hi, position_scale);
    }
    return representable;
}

}

#endif
