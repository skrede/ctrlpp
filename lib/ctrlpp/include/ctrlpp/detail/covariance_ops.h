#ifndef HPP_GUARD_CTRLPP_DETAIL_COVARIANCE_OPS_H
#define HPP_GUARD_CTRLPP_DETAIL_COVARIANCE_OPS_H

/// @brief Covariance-matrix utilities (symmetrisation).
///
/// @cite simon2006 -- Simon, "Optimal State Estimation", 2006, Ch. 5 (covariance conditioning / Joseph form)

#include <Eigen/Core>

namespace ctrlpp::detail
{

/// @brief Symmetrize a square matrix: 0.5 * (M + M^T)
///
/// Numerical operations on covariance matrices can introduce small asymmetries.
/// This enforces exact symmetry, which downstream Cholesky/LDLT decompositions require.
///
/// @cite simon2006 -- Simon, "Optimal State Estimation", 2006, Ch. 5 (numerical covariance conditioning)
template <typename Derived>
[[nodiscard]] inline auto symmetrize(const Eigen::MatrixBase<Derived>& M)
{
    using Scalar = typename Derived::Scalar;
    return (Scalar{0.5} * (M + M.transpose())).eval();
}

}

#endif
