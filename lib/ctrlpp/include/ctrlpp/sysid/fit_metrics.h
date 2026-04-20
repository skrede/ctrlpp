#ifndef HPP_GUARD_CTRLPP_SYSID_FIT_METRICS_H
#define HPP_GUARD_CTRLPP_SYSID_FIT_METRICS_H

/// @brief Model fit metrics (NRMSE, VAF) for system identification validation.
///
/// @cite ljung1999 -- Ljung, "System Identification: Theory for the User", 2nd ed., 1999, Ch. 16 (Model validation)
/// @cite vanoverscheedemoor1996 -- Van Overschee & De Moor, "Subspace Identification for Linear Systems", 1996 (VAF in subspace ID literature)

#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <cstddef>

namespace ctrlpp
{

/// @brief Aggregate of model fit metrics: normalised RMSE and variance accounted for.
///
/// @cite ljung1999 -- Ljung, "System Identification: Theory for the User", 2nd ed., 1999, Ch. 16
template <typename Scalar>
struct fit_metrics
{
    Scalar nrmse{};
    Scalar vaf{};
};

/// @brief Compute NRMSE and VAF for predicted vs actual output trajectories.
///
/// NRMSE is the error norm divided by the centered-output norm (constant-signal fallback to 0 or infinity).
/// VAF is 100 * (1 - var(error) / var(y)) as a percentage (constant-signal fallback to 100 or -infinity).
///
/// @cite ljung1999 -- Ljung, "System Identification: Theory for the User", 2nd ed., 1999, Ch. 16 (Model validation)
/// @cite vanoverscheedemoor1996 -- Van Overschee & De Moor, "Subspace Identification for Linear Systems", 1996 (VAF definition)
template <typename DerivedA, typename DerivedB>
fit_metrics<typename DerivedA::Scalar> compute_fit_metrics(const Eigen::MatrixBase<DerivedA>& y_actual, const Eigen::MatrixBase<DerivedB>& y_predicted)
{
    using Scalar = typename DerivedA::Scalar;

    auto error = (y_actual - y_predicted).eval();
    Scalar mean_y = y_actual.mean();
    auto y_centered = (y_actual.array() - mean_y).matrix().eval();

    Scalar norm_error = error.norm();
    Scalar norm_centered = y_centered.norm();

    // NRMSE
    Scalar nrmse{};
    if(norm_centered < std::numeric_limits<Scalar>::epsilon())
    {
        // Constant y: if prediction is also perfect, NRMSE=0; otherwise infinity
        if(norm_error < std::numeric_limits<Scalar>::epsilon())
            nrmse = Scalar{0};
        else
            nrmse = std::numeric_limits<Scalar>::infinity();
    }
    else
        nrmse = norm_error / norm_centered;

    // VAF
    Scalar var_error{};
    Scalar var_y{};
    auto n = static_cast<Scalar>(y_actual.size());

    if(n > Scalar{1})
    {
        Scalar mean_error = error.mean();
        var_error = (error.array() - mean_error).square().sum() / (n - Scalar{1});
        var_y = (y_actual.array() - mean_y).square().sum() / (n - Scalar{1});
    }

    Scalar vaf{};
    if(var_y < std::numeric_limits<Scalar>::epsilon())
    {
        if(var_error < std::numeric_limits<Scalar>::epsilon())
            vaf = Scalar{100};
        else
            vaf = -std::numeric_limits<Scalar>::infinity();
    }
    else
        vaf = (Scalar{1} - var_error / var_y) * Scalar{100};

    return {.nrmse = nrmse, .vaf = vaf};
}

}

#endif
