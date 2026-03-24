#ifndef HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_JULIER_SIGMA_POINTS_H
#define HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_JULIER_SIGMA_POINTS_H

/// @brief Symmetric sigma point generation (Julier-Uhlmann variant).
///
/// @cite julier2004 -- Julier & Uhlmann, "Unscented Filtering and Nonlinear Estimation", 2004

#include "ctrlpp/types.h"

#include "ctrlpp/estimation/sigma_points/sigma_point_strategy.h"

#include <Eigen/Cholesky>

#include <cmath>
#include <cstddef>

namespace ctrlpp
{

template <typename Scalar>
struct julier_options
{
    Scalar kappa{Scalar{0}};
};

template <typename Scalar, std::size_t NX>
class julier_sigma_points
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr Scalar n = static_cast<Scalar>(NX);

public:
    static constexpr std::size_t num_points = 2 * NX + 1;
    using options_t = julier_options<Scalar>;

    explicit julier_sigma_points(options_t opts = options_t{}) : m_kappa{opts.kappa} {}

    sigma_result<Scalar, NX, num_points> generate(const Vector<Scalar, NX>& x, const Matrix<Scalar, NX, NX>& P) const
    {
        sigma_result<Scalar, NX, num_points> result;

        Scalar gamma = std::sqrt(n + m_kappa);

        // LDLT decomposition of P for matrix square root factor
        Eigen::LDLT<Eigen::Matrix<Scalar, nx, nx>> ldlt(P);

        Eigen::Matrix<Scalar, nx, nx> S;
        if(ldlt.info() == Eigen::Success && ldlt.isPositive())
        {
            // S = L * sqrt(D) where P = L*D*L^T
            Eigen::Matrix<Scalar, nx, nx> L = ldlt.matrixL();
            Vector<Scalar, NX> sqrtD = ldlt.vectorD().cwiseMax(Scalar{0}).cwiseSqrt();
            S = L * sqrtD.asDiagonal();
        }
        else
            S = Eigen::Matrix<Scalar, nx, nx>::Identity();

        // Center sigma point
        result.points[0] = x;

        // Offset sigma points
        for(std::size_t i = 0; i < NX; ++i)
        {
            Vector<Scalar, NX> offset = gamma * S.col(static_cast<int>(i));
            result.points[1 + i] = x + offset;
            result.points[1 + NX + i] = x - offset;
        }

        // Weights: Julier uses Wm = Wc
        Scalar denom = n + m_kappa;
        result.Wm[0] = m_kappa / denom;
        result.Wc[0] = m_kappa / denom;

        Scalar wi = Scalar{1} / (Scalar{2} * denom);
        for(std::size_t i = 1; i < num_points; ++i)
        {
            result.Wm[i] = wi;
            result.Wc[i] = wi;
        }

        return result;
    }

private:
    Scalar m_kappa;
};

} // namespace ctrlpp

#endif
