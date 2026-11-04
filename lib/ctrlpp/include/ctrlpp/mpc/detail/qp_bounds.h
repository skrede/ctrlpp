#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_QP_BOUNDS_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_QP_BOUNDS_H

#include "ctrlpp/mpc/terminal_set.h"

#include "ctrlpp/types.h"

#include "ctrlpp/detail/covariance_ops.h"

#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <variant>
#include <optional>
#include <type_traits>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX>
inline void set_dynamics_bounds(Eigen::VectorX<Scalar>& l, Eigen::VectorX<Scalar>& u, Eigen::Index& row, std::int32_t N, const Vector<Scalar, NX>& x0)
{
    constexpr int nx = static_cast<int>(NX);
    l.segment(row, nx) = x0;
    u.segment(row, nx) = x0;
    row += nx;
    for(std::int32_t k = 1; k <= N; ++k)
    {
        l.segment(row, nx).setZero();
        u.segment(row, nx).setZero();
        row += nx;
    }
}

template <typename Scalar, std::size_t NX>
inline void set_state_bounds(Eigen::VectorX<Scalar>& l, Eigen::VectorX<Scalar>& u, Eigen::Index& row, std::int32_t N, const std::optional<Vector<Scalar, NX>>& x_min,
        const std::optional<Vector<Scalar, NX>>& x_max)
{
    constexpr int nx = static_cast<int>(NX);
    constexpr auto inf = std::numeric_limits<Scalar>::infinity();
    Vector<Scalar, NX> lb = x_min.has_value() ? *x_min : Vector<Scalar, NX>::Constant(-inf);
    Vector<Scalar, NX> ub = x_max.has_value() ? *x_max : Vector<Scalar, NX>::Constant(inf);
    for(std::int32_t k = 0; k < N; ++k)
    {
        l.segment(row, nx) = lb;
        u.segment(row, nx) = ub;
        row += nx;
    }
}

template <typename Scalar, std::size_t NU>
inline void set_input_bounds(Eigen::VectorX<Scalar>& l, Eigen::VectorX<Scalar>& u, Eigen::Index& row, std::int32_t N, const std::optional<Vector<Scalar, NU>>& u_min,
        const std::optional<Vector<Scalar, NU>>& u_max)
{
    constexpr int nu = static_cast<int>(NU);
    constexpr auto inf = std::numeric_limits<Scalar>::infinity();
    Vector<Scalar, NU> lb = u_min.has_value() ? *u_min : Vector<Scalar, NU>::Constant(-inf);
    Vector<Scalar, NU> ub = u_max.has_value() ? *u_max : Vector<Scalar, NU>::Constant(inf);
    for(std::int32_t k = 0; k < N; ++k)
    {
        l.segment(row, nu) = lb;
        u.segment(row, nu) = ub;
        row += nu;
    }
}

template <typename Scalar, std::size_t NU>
inline void set_rate_bounds(Eigen::VectorX<Scalar>& l, Eigen::VectorX<Scalar>& u, Eigen::Index& row, std::int32_t N, const Vector<Scalar, NU>& du_max)
{
    constexpr int nu = static_cast<int>(NU);
    for(std::int32_t k = 0; k < N; ++k)
    {
        l.segment(row, nu) = -du_max;
        u.segment(row, nu) = du_max;
        row += nu;
    }
}

template <typename Scalar, std::size_t NX>
inline void set_ellipsoidal_bounds(Eigen::VectorX<Scalar>& l, Eigen::VectorX<Scalar>& u, Eigen::Index& row, const ellipsoidal_set<Scalar, NX>& s)
{
    constexpr int nx = static_cast<int>(NX);
    const Scalar floor = rounding_bound(nx, s.P.cwiseAbs().maxCoeff(), Scalar{1});
    Eigen::SelfAdjointEigenSolver<Matrix<Scalar, NX, NX>> eig(s.P);
    auto D = eig.eigenvalues();
    // Inscribed box of {x : x^T P x <= alpha} in P's eigenbasis y = V^T x: |y_i| <= sqrt(alpha / (NX lambda_i)) keeps
    // sum_i lambda_i y_i^2 <= alpha. An eigenvalue within the eigensolver's rounding floor of zero (Golub and Van Loan,
    // 4th ed., Sec. 8.1) is a null direction of P, along which the set is open for every alpha >= 0.
    for(Eigen::Index i = 0; i < nx; ++i)
    {
        Scalar bound = D(i) <= floor ? std::numeric_limits<Scalar>::infinity() : std::sqrt(s.alpha / (Scalar{nx} * D(i)));
        l(row + i) = -bound;
        u(row + i) = bound;
    }
    row += nx;
}

template <typename Scalar, std::size_t NX>
inline void set_polytopic_bounds(Eigen::VectorX<Scalar>& l, Eigen::VectorX<Scalar>& u, Eigen::Index& row, const polytopic_set<Scalar, NX>& s)
{
    constexpr Scalar large_bound = Scalar{1e20};
    Eigen::Index n_faces = s.H.rows();
    l.segment(row, n_faces).setConstant(-large_bound);
    u.segment(row, n_faces) = s.h;
    row += n_faces;
}

template <typename Scalar, std::size_t NX>
inline void set_terminal_bounds(Eigen::VectorX<Scalar>& l, Eigen::VectorX<Scalar>& u, Eigen::Index& row, const terminal_set<Scalar, NX>& tset)
{
    std::visit(
            [&](const auto& s)
            {
                using T = std::decay_t<decltype(s)>;
                if constexpr(std::is_same_v<T, ellipsoidal_set<Scalar, NX>>)
                    set_ellipsoidal_bounds<Scalar, NX>(l, u, row, s);
                else
                    set_polytopic_bounds<Scalar, NX>(l, u, row, s);
            },
            tset);
}

}

#endif
