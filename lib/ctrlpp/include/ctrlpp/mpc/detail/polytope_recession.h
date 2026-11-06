#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_POLYTOPE_RECESSION_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_POLYTOPE_RECESSION_H

// A polytope {x : H x <= h} recedes when some nonzero d has H d <= 0. A rank-deficient H leaves a line. Otherwise the
// cone {d : H d <= 0} is pointed, and a pointed cone other than {0} has an extreme ray: the null direction of NX - 1
// independent rows of H, which no row of H has on its positive side. Each row is compared with a candidate direction
// within the NX-term inner-product bound (@cite higham2002, Sec. 3.1) plus the row applied to the direction's forward
// error, the residual of its defining rows mapped back through their pseudo-inverse (as for a vertex, @cite higham2002,
// Sec. 7.7, Eq. 7.31).

#include "ctrlpp/mpc/detail/polytope_redundancy.h"

#include "ctrlpp/mpc/terminal_set.h"

#include "ctrlpp/detail/covariance_ops.h"

#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <optional>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX>
struct recession_candidate
{
    Vector<Scalar, NX> direction;
    Vector<Scalar, NX> error;
};

// The tolerance can only report a direction exact arithmetic might lack, never hide one it has: it errs toward
// refusing a bounded set, which is the sound side.
template <typename Scalar, std::size_t NX>
bool recedes_along(const polytopic_set<Scalar, NX>& set, const recession_candidate<Scalar, NX>& candidate, Scalar sign, Scalar factor)
{
    for(Eigen::Index j = 0; j < set.H.rows(); ++j)
    {
        auto const row = set.H.row(j);
        Scalar const allowance = rounding_bound<Scalar>(static_cast<std::int32_t>(NX), row.cwiseAbs().dot(candidate.direction.cwiseAbs()), factor) + row.cwiseAbs().dot(candidate.error);
        if(sign * row.dot(candidate.direction) > allowance)
            return false;
    }
    return true;
}

// The unit null direction of the rows `subset`, from the last column of Q in the QR of their transpose S^T = Q R. A
// subset whose sigma_min is within the QR's backward error (@cite higham2002, Sec. 19.3) is singular as far as the
// factorization can tell, and determines no direction.
template <typename Scalar, std::size_t NX>
std::optional<recession_candidate<Scalar, NX>> null_direction(const polytopic_set<Scalar, NX>& set, const std::array<Eigen::Index, NX - 1>& subset, Scalar factor)
{
    constexpr int nx = static_cast<int>(NX);
    Eigen::Matrix<Scalar, nx - 1, nx> S;
    for(std::size_t k = 0; k + 1 < NX; ++k)
        S.row(static_cast<Eigen::Index>(k)) = set.H.row(subset[k]);
    Scalar const sigma_min = Eigen::JacobiSVD<Eigen::Matrix<Scalar, nx - 1, nx>>(S).singularValues()(nx - 2);
    if(sigma_min <= rounding_bound<Scalar>(static_cast<std::int32_t>(NX * (NX - 1)), S.norm(), factor))
        return std::nullopt;
    Eigen::HouseholderQR<Eigen::Matrix<Scalar, nx, nx - 1>> const qr(S.transpose());
    Matrix<Scalar, NX, NX> const Q = qr.householderQ();
    Eigen::Matrix<Scalar, nx - 1, nx - 1> const R = qr.matrixQR().template topRows<nx - 1>().template triangularView<Eigen::Upper>();
    Eigen::Matrix<Scalar, nx - 1, nx - 1> const R_inv_T = R.transpose().template triangularView<Eigen::Lower>().solve(Eigen::Matrix<Scalar, nx - 1, nx - 1>::Identity());
    recession_candidate<Scalar, NX> candidate{.direction = Q.col(nx - 1), .error = Vector<Scalar, NX>::Zero()};
    Eigen::Matrix<Scalar, nx - 1, 1> slack;
    for(Eigen::Index k = 0; k < nx - 1; ++k)
        slack(k) = std::abs(S.row(k).dot(candidate.direction)) + rounding_bound<Scalar>(static_cast<std::int32_t>(NX), S.row(k).cwiseAbs().dot(candidate.direction.cwiseAbs()), factor);
    candidate.error = (Q.template leftCols<nx - 1>() * R_inv_T).cwiseAbs() * slack;
    return candidate;
}

template <typename Scalar, std::size_t NX>
bool has_extreme_ray(const polytopic_set<Scalar, NX>& set, Scalar factor)
{
    std::array<Eigen::Index, NX - 1> subset{};
    std::iota(subset.begin(), subset.end(), Eigen::Index{0});
    do
    {
        auto const candidate = null_direction<Scalar, NX>(set, subset, factor);
        if(candidate && (recedes_along<Scalar, NX>(set, *candidate, Scalar{1}, factor) || recedes_along<Scalar, NX>(set, *candidate, Scalar{-1}, factor)))
            return true;
    } while(next_subset<NX - 1>(subset, set.H.rows()));
    return false;
}

template <typename Scalar, std::size_t NX>
bool has_recession_direction(const polytopic_set<Scalar, NX>& set, Scalar factor)
{
    constexpr auto nx = static_cast<Eigen::Index>(NX);
    if(set.H.rows() < nx || set.H.colPivHouseholderQr().rank() < nx)
        return true;
    if constexpr(NX == 1)
    {
        recession_candidate<Scalar, NX> const unit{.direction = Vector<Scalar, NX>::Ones(), .error = Vector<Scalar, NX>::Zero()};
        return recedes_along<Scalar, NX>(set, unit, Scalar{1}, factor) || recedes_along<Scalar, NX>(set, unit, Scalar{-1}, factor);
    }
    else
        return has_extreme_ray<Scalar, NX>(set, factor);
}

}

#endif
