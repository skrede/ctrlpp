#ifndef HPP_GUARD_CTRLPP_CONTROL_DETAIL_LQR_CONTINUOUS_H
#define HPP_GUARD_CTRLPP_CONTROL_DETAIL_LQR_CONTINUOUS_H

#include "ctrlpp/control/care_types.h"

#include "ctrlpp/expected.h"

#include <Eigen/Dense>

#include <limits>
#include <cstddef>

namespace ctrlpp::detail
{

// Eigen's LDLT solve zeroes rank-deficient directions instead of producing infinities, so a singular R would pass
// every finiteness check; the pivots take the Riccati builds' reciprocal-pivot test, and sign is not tested. An R just
// above that threshold gives a large R^{-1}, and B R^{-1} B^T can still overflow.
template <typename Scalar, std::size_t NX, std::size_t NU>
ctrlpp::expected<Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>, care_error> lqr_continuous_hamiltonian(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
        const Eigen::Matrix<Scalar, int(NX), int(NU)>& B, const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q_scaled, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R_scaled,
        Eigen::Matrix<Scalar, int(NU), int(NU)>& R_inv)
{
    constexpr int nx = static_cast<int>(NX);

    const auto R_ldlt = R_scaled.ldlt();
    const auto R_pivots = R_ldlt.vectorD().cwiseAbs().eval();
    if(!(R_pivots.minCoeff() > Scalar{static_cast<int>(NU)} * std::numeric_limits<Scalar>::epsilon() * R_pivots.maxCoeff()))
        return ctrlpp::unexpected(care_error::singular_r);
    R_inv = R_ldlt.solve(Eigen::Matrix<Scalar, int(NU), int(NU)>::Identity());

    Eigen::Matrix<Scalar, 2 * nx, 2 * nx> H;
    H.template block<nx, nx>(0, 0) = A;
    H.template block<nx, nx>(0, nx) = -B * R_inv * B.transpose();
    H.template block<nx, nx>(nx, 0) = -Q_scaled;
    H.template block<nx, nx>(nx, nx) = -A.transpose();
    if(!H.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);
    return H;
}

}

#endif
