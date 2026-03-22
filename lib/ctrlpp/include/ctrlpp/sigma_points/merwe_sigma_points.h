#ifndef HPP_GUARD_CTRLPP_SIGMA_POINTS_MERWE_SIGMA_POINTS_H
#define HPP_GUARD_CTRLPP_SIGMA_POINTS_MERWE_SIGMA_POINTS_H

#include "ctrlpp/sigma_points/sigma_point_strategy.h"
#include "ctrlpp/types.h"

#include <Eigen/Cholesky>

#include <cmath>
#include <cstddef>

namespace ctrlpp {

template<typename Scalar>
struct merwe_options {
    Scalar alpha{Scalar{1e-3}};
    Scalar beta{Scalar{2}};
    Scalar kappa{Scalar{0}};
};

template<typename Scalar, std::size_t NX>
class merwe_sigma_points {
    static constexpr int nx = static_cast<int>(NX);
    static constexpr Scalar n = static_cast<Scalar>(NX);

public:
    static constexpr std::size_t num_points = 2 * NX + 1;
    using options_t = merwe_options<Scalar>;

    explicit merwe_sigma_points(options_t opts = options_t{})
        : alpha_{opts.alpha}
        , beta_{opts.beta}
        , kappa_{opts.kappa}
    {
    }

    [[nodiscard]] auto generate(const Vector<Scalar, NX>& x,
                                const Matrix<Scalar, NX, NX>& P) const
        -> sigma_result<Scalar, NX, num_points>
    {
        sigma_result<Scalar, NX, num_points> result;

        Scalar lambda = alpha_ * alpha_ * (n + kappa_) - n;
        Scalar gamma = std::sqrt(n + lambda);

        // LDLT decomposition of P for matrix square root factor
        Eigen::LDLT<Eigen::Matrix<Scalar, nx, nx>> ldlt(P);

        Eigen::Matrix<Scalar, nx, nx> S;
        if (ldlt.info() == Eigen::Success && ldlt.isPositive()) {
            // S = L * sqrt(D) where P = L*D*L^T
            Eigen::Matrix<Scalar, nx, nx> L = ldlt.matrixL();
            Vector<Scalar, NX> sqrtD = ldlt.vectorD().cwiseMax(Scalar{0}).cwiseSqrt();
            S = L * sqrtD.asDiagonal();
        } else {
            // Fallback: identity scaling (defensive)
            S = Eigen::Matrix<Scalar, nx, nx>::Identity();
        }

        // Center sigma point
        result.points[0] = x;

        // Offset sigma points: x +/- gamma * columns of S
        for (std::size_t i = 0; i < NX; ++i) {
            Vector<Scalar, NX> offset = gamma * S.col(static_cast<int>(i));
            result.points[1 + i] = x + offset;
            result.points[1 + NX + i] = x - offset;
        }

        // Weights
        Scalar denom = n + lambda;
        result.Wm[0] = lambda / denom;
        result.Wc[0] = lambda / denom + (Scalar{1} - alpha_ * alpha_ + beta_);

        Scalar wi = Scalar{1} / (Scalar{2} * denom);
        for (std::size_t i = 1; i < num_points; ++i) {
            result.Wm[i] = wi;
            result.Wc[i] = wi;
        }

        return result;
    }

private:
    Scalar alpha_;
    Scalar beta_;
    Scalar kappa_;
};

}

#endif
