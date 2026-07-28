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

/// @brief Resolution floor below which a norm formed from an n-sample record
/// carries no significant digits.
///
/// Both fallbacks below ask the same question -- has this vector vanished? -- and
/// the answer depends on the record's own magnitude, never on an absolute
/// constant. An absolute one is wrong in both directions and both are reachable:
/// at a record scale of 1e-17 it calls a genuinely varying record constant and
/// reports a perfect fit for a predictor that explains nothing, and at a large
/// scale it calls a constant record varying and then divides one rounding-level
/// quantity by another.
///
/// The count is enumerated rather than chosen. Forming the mean costs n - 1
/// additions and one division, that is n operations, each worth one unit in the
/// last place at the scale of the largest sample; each centered entry carries
/// those plus its own subtraction, so n + 1; and the Euclidean norm of n entries
/// each bounded by that is at most the square root of n times it. Every operation
/// is counted whether or not it rounds, so the floor bounds the accumulated
/// departure from above rather than describing it tightly.
template <typename Derived>
auto norm_resolution_floor(const Eigen::MatrixBase<Derived>& y) -> typename Derived::Scalar
{
    using Scalar = typename Derived::Scalar;
    auto const n = static_cast<Scalar>(y.size());
    Scalar const largest = y.cwiseAbs().maxCoeff();
    return std::sqrt(n) * (n + Scalar{1}) * std::numeric_limits<Scalar>::epsilon() * largest;
}

/// @brief Compute NRMSE and VAF for predicted vs actual output trajectories.
///
/// NRMSE is the error norm divided by the centered-output norm (constant-signal fallback to 0 or infinity).
/// VAF is 100 * (1 - var(error) / var(y)) as a percentage (constant-signal fallback to 100 or -infinity).
///
/// Whether either quantity counts as vanishing is decided against
/// `norm_resolution_floor` -- a counted rounding chain at the scale of the record
/// that entered -- and never against an absolute constant.
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

    Scalar const norm_floor = norm_resolution_floor(y_actual);

    // NRMSE
    Scalar nrmse{};
    if(norm_centered <= norm_floor)
    {
        // Constant y: if prediction is also perfect, NRMSE=0; otherwise infinity
        if(norm_error <= norm_floor)
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

    // The same floor expressed in the units the variances carry: a norm squared
    // over the same degrees of freedom the variances are divided by.
    Scalar const var_floor = (n > Scalar{1}) ? norm_floor * norm_floor / (n - Scalar{1})
                                             : norm_floor * norm_floor;

    Scalar vaf{};
    if(var_y <= var_floor)
    {
        if(var_error <= var_floor)
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
