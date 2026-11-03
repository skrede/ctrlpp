#ifndef HPP_GUARD_CTRLPP_SYSID_FIT_METRICS_H
#define HPP_GUARD_CTRLPP_SYSID_FIT_METRICS_H

// @cite ljung1999 -- Ljung, "System Identification: Theory for the User", 2nd ed., 1999, Ch. 16 (Model validation)
// @cite vanoverscheedemoor1996 -- Van Overschee & De Moor, "Subspace Identification for Linear Systems", 1996 (VAF in subspace ID literature)

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <cstddef>

namespace ctrlpp
{

template <typename Scalar>
struct fit_metrics
{
    fit_metrics()
        : nrmse{}
        , vaf{}
    {
    }

    Scalar nrmse;
    Scalar vaf;
};

enum class fit_metrics_error
{
    empty_record,
    size_mismatch,
    non_finite_record,
    not_column_vector,
};

// The floor counts every operation of the mean, the centering and the norm, sqrt(n) * (n + 1) units in the last
// place at the scale of the largest sample, so whether a norm has vanished is judged against the record's own magnitude.
template <typename Derived>
typename Derived::Scalar norm_resolution_floor(const Eigen::MatrixBase<Derived>& y)
{
    using Scalar = typename Derived::Scalar;
    if(y.size() == 0)
        return Scalar{0};
    auto const n = static_cast<Scalar>(y.size());
    Scalar const largest = y.cwiseAbs().maxCoeff();
    return std::sqrt(n) * (n + Scalar{1}) * std::numeric_limits<Scalar>::epsilon() * largest;
}

namespace detail
{

template <typename Scalar>
Scalar normalized_rmse(Scalar norm_error, Scalar norm_centered, Scalar norm_floor)
{
    if(norm_centered <= norm_floor)
        return norm_error <= norm_floor ? Scalar{0} : std::numeric_limits<Scalar>::infinity();
    return norm_error / norm_centered;
}

template <typename Scalar>
Scalar variance_accounted_for(Scalar var_error, Scalar var_y, Scalar norm_floor, Scalar n)
{
    // The norm floor in the units the variances carry: a norm squared over the variances' degrees of freedom.
    Scalar const var_floor = (n > Scalar{1}) ? norm_floor * norm_floor / (n - Scalar{1}) : norm_floor * norm_floor;
    if(var_y <= var_floor)
        return var_error <= var_floor ? Scalar{100} : -std::numeric_limits<Scalar>::infinity();
    return (Scalar{1} - var_error / var_y) * Scalar{100};
}

template <typename DerivedA, typename DerivedB>
fit_metrics<typename DerivedA::Scalar> compute_fit_metrics_unchecked(const Eigen::MatrixBase<DerivedA>& y_actual, const Eigen::MatrixBase<DerivedB>& y_predicted)
{
    using Scalar = typename DerivedA::Scalar;

    auto error = (y_actual - y_predicted).eval();
    Scalar mean_y = y_actual.mean();
    auto y_centered = (y_actual.array() - mean_y).matrix().eval();
    Scalar norm_error = error.norm();
    Scalar norm_centered = y_centered.norm();
    Scalar const norm_floor = norm_resolution_floor(y_actual);

    fit_metrics<Scalar> metrics;
    metrics.nrmse = normalized_rmse(norm_error, norm_centered, norm_floor);

    Scalar var_error{};
    Scalar var_y{};
    auto n = static_cast<Scalar>(y_actual.size());
    if(n > Scalar{1})
    {
        Scalar mean_error = error.mean();
        var_error = (error.array() - mean_error).square().sum() / (n - Scalar{1});
        var_y = (y_actual.array() - mean_y).square().sum() / (n - Scalar{1});
    }
    metrics.vaf = variance_accounted_for(var_error, var_y, norm_floor, n);
    return metrics;
}

}

template <typename DerivedA, typename DerivedB>
ctrlpp::expected<fit_metrics<typename DerivedA::Scalar>, fit_metrics_error> compute_fit_metrics(
        const Eigen::MatrixBase<DerivedA>& y_actual, const Eigen::MatrixBase<DerivedB>& y_predicted)
{
    if(y_actual.size() == 0 || y_predicted.size() == 0)
        return ctrlpp::unexpected(fit_metrics_error::empty_record);
    if(y_actual.rows() != y_predicted.rows() || y_actual.cols() != y_predicted.cols())
        return ctrlpp::unexpected(fit_metrics_error::size_mismatch);
    if(y_actual.cols() != 1 || y_predicted.cols() != 1)
        return ctrlpp::unexpected(fit_metrics_error::not_column_vector);
    if(!y_actual.allFinite() || !y_predicted.allFinite())
        return ctrlpp::unexpected(fit_metrics_error::non_finite_record);
    return detail::compute_fit_metrics_unchecked(y_actual, y_predicted);
}

}

#endif
