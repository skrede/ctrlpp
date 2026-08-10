#ifndef HPP_GUARD_CTRLPP_SYSID_SYSID_RESULT_H
#define HPP_GUARD_CTRLPP_SYSID_SYSID_RESULT_H

#include "ctrlpp/model/state_space.h"

#include "ctrlpp/sysid/fit_metrics.h"

#include <Eigen/Dense>

namespace ctrlpp
{

/// @brief What the least-squares stage of an ARX identification resolved, as
/// distinct from how well the identified model reproduces the record.
template <typename Scalar>
struct arx_diagnostics
{
    /// Regressor rows the fit was formed from: the record length less
    /// max(NA, NB). The goodness-of-fit metrics score every sample of the
    /// record instead, so the two counts differ by the startup transient.
    std::size_t effective_samples{};
    /// Regressor columns, NA + NB.
    std::size_t parameter_count{};
    /// Rank of the regressor matrix, counted by the column-pivoting QR against
    /// the active threshold. Below `parameter_count` the fit is rank-deficient
    /// and the coefficients in the unresolved directions are not determined by
    /// the record.
    std::size_t numerical_rank{};
    /// Euclidean norm of the least-squares residual over the fitted rows. This
    /// is the residual of the FIT, over the same rows the rank was counted
    /// from, not the residual of the simulation the metrics score.
    Scalar residual_norm{};
};

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct arx_result
{
    discrete_state_space<Scalar, NX, NU, NY> system{};
    fit_metrics<Scalar> metrics{};
    arx_diagnostics<Scalar> diagnostics{};
};

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct moesp_result
{
    discrete_state_space<Scalar, NX, NU, NY> system{};
    Eigen::VectorX<Scalar> singular_values{};
    fit_metrics<Scalar> metrics{};
    Scalar condition_number{};
};

}

#endif
