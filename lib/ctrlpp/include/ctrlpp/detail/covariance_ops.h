#ifndef HPP_GUARD_CTRLPP_DETAIL_COVARIANCE_OPS_H
#define HPP_GUARD_CTRLPP_DETAIL_COVARIANCE_OPS_H

#include <Eigen/Core>

namespace ctrlpp::detail
{

/// @brief Symmetrize a square matrix: 0.5 * (M + M^T)
///
/// Numerical operations on covariance matrices can introduce small asymmetries.
/// This enforces exact symmetry, which downstream Cholesky/LDLT decompositions require.
template <typename Derived>
[[nodiscard]] inline auto symmetrize(const Eigen::MatrixBase<Derived>& M)
{
    using Scalar = typename Derived::Scalar;
    return (Scalar{0.5} * (M + M.transpose())).eval();
}

}

#endif
