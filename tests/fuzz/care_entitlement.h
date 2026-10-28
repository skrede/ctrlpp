#ifndef HPP_GUARD_CTRLPP_TESTS_FUZZ_CARE_ENTITLEMENT_H
#define HPP_GUARD_CTRLPP_TESTS_FUZZ_CARE_ENTITLEMENT_H

#include "care_quad_reference.h"

#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <optional>

namespace ctrlpp::fuzz
{

namespace detail
{

inline auto kron(const Eigen::Matrix2d &left, const Eigen::Matrix2d &right) -> Eigen::Matrix4d
{
    Eigen::Matrix4d product;
    for(Eigen::Index i = 0; i < 2; ++i)
        for(Eigen::Index j = 0; j < 2; ++j)
            product.block<2, 2>(2 * i, 2 * j) = left(i, j) * right;
    return product;
}

// vec(W') = K vec(W) under column-major vec, which swaps the two off-diagonal entries.
inline auto commutation() -> Eigen::Matrix4d
{
    Eigen::Matrix4d K = Eigen::Matrix4d::Zero();
    K(0, 0)           = 1.0;
    K(1, 2)           = 1.0;
    K(2, 1)           = 1.0;
    K(3, 3)           = 1.0;
    return K;
}

inline auto to_matrix(const double M[2][2]) -> Eigen::Matrix2d
{
    Eigen::Matrix2d out;
    out << M[0][0], M[0][1], M[1][0], M[1][1];
    return out;
}

inline auto to_matrix(const quad M[2][2]) -> Eigen::Matrix2d
{
    const double rounded[2][2] = {{static_cast<double>(M[0][0]), static_cast<double>(M[0][1])}, {static_cast<double>(M[1][0]), static_cast<double>(M[1][1])}};
    return to_matrix(rounded);
}

// JacobiSVD leaves its singular values unset on a non-finite input and reports it only in info().
inline auto largest_singular_value(const Eigen::Matrix4d &M) -> double
{
    const Eigen::JacobiSVD<Eigen::Matrix4d> svd(M);
    return svd.info() == Eigen::Success ? svd.singularValues()(0) : std::numeric_limits<double>::quiet_NaN();
}

}

/// @brief The condition number of the continuous Riccati equation A'X + XA + Q - XGX = 0,
/// G = B R^-1 B', at its stabilizing solution X:
/// `(||Theta|| ||A|| + ||Q|| / sep + ||Pi|| ||G||) / ||X||`, every operator norm induced by the
/// Frobenius norm and evaluated exactly on its 4x4 Kronecker matrix.
///
/// X is the converged binary128 reference rounded to binary64; nothing else of the solve is read.
/// Returns no value when the reference has not converged, X is zero, the closed loop's Lyapunov
/// operator is singular, or any quantity is not finite.
///
/// @cite petkov1999 -- Petkov and Sima, "SB02QD", SLICOT Library, 1999, METHOD
inline auto care_condition_number(const double A[2][2], const double B[2], const double Q[2][2], double R, const quad_reference_result &reference) -> std::optional<double>
{
    if(!reference.converged)
        return std::nullopt;
    const Eigen::Matrix2d identity = Eigen::Matrix2d::Identity();
    const Eigen::Matrix2d a        = detail::to_matrix(A);
    const Eigen::Matrix2d x        = detail::to_matrix(reference.P);
    const Eigen::Vector2d b{B[0], B[1]};
    const Eigen::Matrix2d g             = b * b.transpose() / R;
    const Eigen::Matrix2d closed_loop_t = (a - g * x).transpose();
    const Eigen::Matrix4d omega         = detail::kron(identity, closed_loop_t) + detail::kron(closed_loop_t, identity);
    if(!omega.allFinite() || !x.allFinite() || !(x.norm() > 0.0))
        return std::nullopt;
    const Eigen::JacobiSVD<Eigen::Matrix4d> omega_svd(omega, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const double separation = omega_svd.singularValues()(3);
    if(!(separation > 0.0))
        return std::nullopt;
    const Eigen::Matrix4d omega_inverse = omega_svd.solve(Eigen::Matrix4d::Identity());
    const Eigen::Matrix4d theta         = omega_inverse * (detail::kron(x, identity) * detail::commutation() + detail::kron(identity, x));
    const Eigen::Matrix4d pi            = omega_inverse * detail::kron(x, x);
    const double condition = (detail::largest_singular_value(theta) * a.norm() + detail::to_matrix(Q).norm() / separation + detail::largest_singular_value(pi) * g.norm()) / x.norm();
    if(!std::isfinite(condition))
        return std::nullopt;
    return condition;
}

/// @brief Whether a backward-stable binary64 answer on this pose is entitled to a relative forward
/// error of at most sqrt(eps), half a significand.
///
/// `backward_error_factor` assumes the answer solves a problem within that many units of relative
/// backward error, so its forward error is at most about that factor times the condition number
/// times eps. The verdict reads only the reference and the posed data.
inline auto care_forward_error_entitled(const double A[2][2], const double B[2], const double Q[2][2], double R, const quad_reference_result &reference,
                                        double backward_error_factor = 1.0) -> bool
{
    constexpr double eps                  = std::numeric_limits<double>::epsilon();
    const std::optional<double> condition = care_condition_number(A, B, Q, R, reference);
    return condition.has_value() && backward_error_factor * *condition * eps <= std::sqrt(eps);
}

}

#endif
