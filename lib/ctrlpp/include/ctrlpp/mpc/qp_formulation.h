#ifndef HPP_GUARD_CTRLPP_MPC_QP_FORMULATION_H
#define HPP_GUARD_CTRLPP_MPC_QP_FORMULATION_H

/// @brief Sparse QP formulation for linear MPC (condensed/sparse form).
///
/// @cite rawlings2017 -- Rawlings et al., "Model Predictive Control: Theory, Computation, and Design", 2017, Ch. 2

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/terminal_set.h"
#include "ctrlpp/types.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace ctrlpp::detail
{

// --- Terminal constraint row count helper ---

template <typename Scalar, std::size_t NX>
inline int terminal_constraint_rows(const std::optional<terminal_set<Scalar, NX>>& terminal_constraint)
{
    constexpr int nx = static_cast<int>(NX);
    if(!terminal_constraint.has_value())
        return 0;
    return std::visit(
        [](const auto& s) -> int
        {
            using T = std::decay_t<decltype(s)>;
            if constexpr(std::is_same_v<T, ellipsoidal_set<Scalar, NX>>)
                return 2 * nx;
            else
                return static_cast<int>(s.H.rows());
        },
        terminal_constraint.value());
}

// --- Cost matrix building ---

template <typename Scalar, std::size_t NX, std::size_t NU>
inline void add_diagonal_block(std::vector<Eigen::Triplet<Scalar>>& trips, int offset, const auto& block, int size)
{
    for(int i = 0; i < size; ++i)
        for(int j = i; j < size; ++j)
            if(block(i, j) != Scalar{0})
                trips.emplace_back(offset + i, offset + j, block(i, j));
}

template <typename Scalar, std::size_t NX, std::size_t NU>
inline void add_soft_constraint_triplets(std::vector<Eigen::Triplet<Scalar>>& trips, int N, int slack_offset,
                                         Scalar soft_penalty, const std::optional<Vector<Scalar, NX>>& soft_state_penalty)
{
    constexpr int nx = static_cast<int>(NX);
    for(int k = 0; k < N; ++k)
        for(int i = 0; i < nx; ++i)
        {
            Scalar penalty = soft_state_penalty.has_value() ? (*soft_state_penalty)(i) : soft_penalty;
            trips.emplace_back(slack_offset + k * nx + i, slack_offset + k * nx + i, penalty);
        }
}

template <typename Scalar, std::size_t NX, std::size_t NU>
[[nodiscard]] auto build_cost_matrix(int N,
                                     const Matrix<Scalar, NX, NX>& Q,
                                     const Matrix<Scalar, NU, NU>& R,
                                     const Matrix<Scalar, NX, NX>& Qf,
                                     bool has_soft_constraints,
                                     Scalar soft_penalty,
                                     const std::optional<Vector<Scalar, NX>>& soft_state_penalty = {}) -> Eigen::SparseMatrix<Scalar, Eigen::ColMajor>
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);

    int n_x_total = (N + 1) * nx;
    int n_u_total = N * nu;
    int n_slack = has_soft_constraints ? N * nx : 0;
    int n_dec = n_x_total + n_u_total + n_slack;

    std::vector<Eigen::Triplet<Scalar>> triplets;
    triplets.reserve(static_cast<std::size_t>((N + 1) * nx * (nx + 1) / 2 + N * nu * (nu + 1) / 2 + n_slack));

    for(int k = 0; k < N; ++k)
        add_diagonal_block<Scalar, NX, NU>(triplets, k * nx, Q, nx);
    add_diagonal_block<Scalar, NX, NU>(triplets, N * nx, Qf, nx);
    for(int k = 0; k < N; ++k)
        add_diagonal_block<Scalar, NX, NU>(triplets, n_x_total + k * nu, R, nu);
    if(has_soft_constraints)
        add_soft_constraint_triplets<Scalar, NX, NU>(triplets, N, n_x_total + n_u_total, soft_penalty, soft_state_penalty);

    Eigen::SparseMatrix<Scalar, Eigen::ColMajor> P(n_dec, n_dec);
    P.setFromTriplets(triplets.begin(), triplets.end());
    return P;
}

// --- Constraint matrix building: sub-steps ---

template <typename Scalar, std::size_t NX, std::size_t NU>
inline void add_dynamics_triplets(std::vector<Eigen::Triplet<Scalar>>& trips,
                                  int& row,
                                  int N,
                                  int n_x_total,
                                  const Matrix<Scalar, NX, NX>& A_sys,
                                  const Matrix<Scalar, NX, NU>& B_sys)
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);

    for(int i = 0; i < nx; ++i)
        trips.emplace_back(row + i, i, Scalar{1});
    row += nx;

    for(int k = 1; k <= N; ++k)
    {
        int x_k = k * nx;
        int x_km1 = (k - 1) * nx;
        int u_km1 = n_x_total + (k - 1) * nu;

        for(int i = 0; i < nx; ++i)
            for(int j = 0; j < nx; ++j)
                if(A_sys(i, j) != Scalar{0})
                    trips.emplace_back(row + i, x_km1 + j, -A_sys(i, j));

        for(int i = 0; i < nx; ++i)
            trips.emplace_back(row + i, x_k + i, Scalar{1});

        for(int i = 0; i < nx; ++i)
            for(int j = 0; j < nu; ++j)
                if(B_sys(i, j) != Scalar{0})
                    trips.emplace_back(row + i, u_km1 + j, -B_sys(i, j));

        row += nx;
    }
}

template <typename Scalar, std::size_t NX, std::size_t NU>
inline void add_state_bound_triplets(std::vector<Eigen::Triplet<Scalar>>& trips, int& row, int N, int n_x_total, int n_u_total, bool has_soft)
{
    constexpr int nx = static_cast<int>(NX);
    for(int k = 1; k <= N; ++k)
    {
        int x_k = k * nx;
        for(int i = 0; i < nx; ++i)
            trips.emplace_back(row + i, x_k + i, Scalar{1});

        if(has_soft)
        {
            int slack = n_x_total + n_u_total + (k - 1) * nx;
            for(int i = 0; i < nx; ++i)
                trips.emplace_back(row + i, slack + i, Scalar{-1});
        }
        row += nx;
    }
}

template <typename Scalar, std::size_t NX, std::size_t NU>
inline void add_input_bound_triplets(std::vector<Eigen::Triplet<Scalar>>& trips, int& row, int N, int n_x_total)
{
    constexpr int nu = static_cast<int>(NU);
    for(int k = 0; k < N; ++k)
    {
        int u_k = n_x_total + k * nu;
        for(int i = 0; i < nu; ++i)
            trips.emplace_back(row + i, u_k + i, Scalar{1});
        row += nu;
    }
}

template <typename Scalar, std::size_t NX, std::size_t NU>
inline void add_rate_bound_triplets(std::vector<Eigen::Triplet<Scalar>>& trips, int& row, int N, int n_x_total)
{
    constexpr int nu = static_cast<int>(NU);
    for(int k = 0; k < N; ++k)
    {
        int u_k = n_x_total + k * nu;
        for(int i = 0; i < nu; ++i)
            trips.emplace_back(row + i, u_k + i, Scalar{1});
        if(k > 0)
        {
            int u_km1 = n_x_total + (k - 1) * nu;
            for(int i = 0; i < nu; ++i)
                trips.emplace_back(row + i, u_km1 + i, Scalar{-1});
        }
        row += nu;
    }
}

template <typename Scalar, std::size_t NX>
inline void add_ellipsoidal_triplets(std::vector<Eigen::Triplet<Scalar>>& trips, int& row, int x_N_offset, const ellipsoidal_set<Scalar, NX>& s)
{
    constexpr int nx = static_cast<int>(NX);
    Eigen::SelfAdjointEigenSolver<Matrix<Scalar, NX, NX>> eig(s.P);
    auto V = eig.eigenvectors();

    for(int i = 0; i < nx; ++i)
        for(int j = 0; j < nx; ++j)
            if(V(j, i) != Scalar{0})
            {
                trips.emplace_back(row + 2 * i, x_N_offset + j, V(j, i));
                trips.emplace_back(row + 2 * i + 1, x_N_offset + j, -V(j, i));
            }
    row += 2 * nx;
}

template <typename Scalar, std::size_t NX>
inline void add_polytopic_triplets(std::vector<Eigen::Triplet<Scalar>>& trips, int& row, int x_N_offset, const polytopic_set<Scalar, NX>& s)
{
    constexpr int nx = static_cast<int>(NX);
    int n_faces = static_cast<int>(s.H.rows());
    for(int i = 0; i < n_faces; ++i)
        for(int j = 0; j < nx; ++j)
            if(s.H(i, j) != Scalar{0})
                trips.emplace_back(row + i, x_N_offset + j, s.H(i, j));
    row += n_faces;
}

template <typename Scalar, std::size_t NX>
inline void add_terminal_set_triplets(std::vector<Eigen::Triplet<Scalar>>& trips, int& row, int x_N_offset, const terminal_set<Scalar, NX>& tset)
{
    std::visit([&](const auto& s) {
        using T = std::decay_t<decltype(s)>;
        if constexpr(std::is_same_v<T, ellipsoidal_set<Scalar, NX>>)
            add_ellipsoidal_triplets<Scalar, NX>(trips, row, x_N_offset, s);
        else
            add_polytopic_triplets<Scalar, NX>(trips, row, x_N_offset, s);
    }, tset);
}

template <typename Scalar, std::size_t NX, std::size_t NU>
[[nodiscard]] auto build_constraint_matrix(int N,
                                           const Matrix<Scalar, NX, NX>& A_sys,
                                           const Matrix<Scalar, NX, NU>& B_sys,
                                           bool has_state_bounds,
                                           bool has_soft_constraints,
                                           bool has_input_bounds,
                                           bool has_rate_bounds,
                                           const std::optional<terminal_set<Scalar, NX>>& terminal_constraint = {}) -> Eigen::SparseMatrix<Scalar, Eigen::ColMajor>
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);

    int n_x_total = (N + 1) * nx;
    int n_u_total = N * nu;
    int n_slack = has_soft_constraints ? N * nx : 0;
    int n_dec = n_x_total + n_u_total + n_slack;
    int n_con = (N + 1) * nx + (has_state_bounds ? N * nx : 0) + (has_input_bounds ? N * nu : 0) + (has_rate_bounds ? N * nu : 0) + terminal_constraint_rows<Scalar, NX>(terminal_constraint);

    std::vector<Eigen::Triplet<Scalar>> triplets;
    triplets.reserve(static_cast<std::size_t>((N + 1) * nx + N * nx * nx + N * nx * nu + (has_state_bounds ? N * nx : 0) + n_slack + (has_input_bounds ? N * nu : 0) + (has_rate_bounds ? 2 * N * nu : 0)));

    int row = 0;
    add_dynamics_triplets<Scalar, NX, NU>(triplets, row, N, n_x_total, A_sys, B_sys);
    if(has_state_bounds)
        add_state_bound_triplets<Scalar, NX, NU>(triplets, row, N, n_x_total, n_u_total, has_soft_constraints);
    if(has_input_bounds)
        add_input_bound_triplets<Scalar, NX, NU>(triplets, row, N, n_x_total);
    if(has_rate_bounds)
        add_rate_bound_triplets<Scalar, NX, NU>(triplets, row, N, n_x_total);
    if(terminal_constraint.has_value())
        add_terminal_set_triplets<Scalar, NX>(triplets, row, N * nx, terminal_constraint.value());

    Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A(n_con, n_dec);
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A;
}

// --- Bounds vector building: sub-steps ---

template <typename Scalar, std::size_t NX>
inline void set_dynamics_bounds(Eigen::VectorX<Scalar>& l, Eigen::VectorX<Scalar>& u, int& row, int N, const Vector<Scalar, NX>& x0)
{
    constexpr int nx = static_cast<int>(NX);
    l.segment(row, nx) = x0;
    u.segment(row, nx) = x0;
    row += nx;
    for(int k = 1; k <= N; ++k)
    {
        l.segment(row, nx).setZero();
        u.segment(row, nx).setZero();
        row += nx;
    }
}

template <typename Scalar, std::size_t NX>
inline void set_state_bounds(Eigen::VectorX<Scalar>& l, Eigen::VectorX<Scalar>& u, int& row, int N,
                             const std::optional<Vector<Scalar, NX>>& x_min, const std::optional<Vector<Scalar, NX>>& x_max)
{
    constexpr int nx = static_cast<int>(NX);
    constexpr auto inf = std::numeric_limits<Scalar>::infinity();
    Vector<Scalar, NX> lb = x_min.has_value() ? x_min.value() : Vector<Scalar, NX>::Constant(-inf);
    Vector<Scalar, NX> ub = x_max.has_value() ? x_max.value() : Vector<Scalar, NX>::Constant(inf);
    for(int k = 0; k < N; ++k)
    {
        l.segment(row, nx) = lb;
        u.segment(row, nx) = ub;
        row += nx;
    }
}

template <typename Scalar, std::size_t NU>
inline void set_input_bounds(Eigen::VectorX<Scalar>& l, Eigen::VectorX<Scalar>& u, int& row, int N,
                             const std::optional<Vector<Scalar, NU>>& u_min, const std::optional<Vector<Scalar, NU>>& u_max)
{
    constexpr int nu = static_cast<int>(NU);
    constexpr auto inf = std::numeric_limits<Scalar>::infinity();
    Vector<Scalar, NU> lb = u_min.has_value() ? u_min.value() : Vector<Scalar, NU>::Constant(-inf);
    Vector<Scalar, NU> ub = u_max.has_value() ? u_max.value() : Vector<Scalar, NU>::Constant(inf);
    for(int k = 0; k < N; ++k)
    {
        l.segment(row, nu) = lb;
        u.segment(row, nu) = ub;
        row += nu;
    }
}

template <typename Scalar, std::size_t NU>
inline void set_rate_bounds(Eigen::VectorX<Scalar>& l, Eigen::VectorX<Scalar>& u, int& row, int N, const Vector<Scalar, NU>& du_max)
{
    constexpr int nu = static_cast<int>(NU);
    for(int k = 0; k < N; ++k)
    {
        l.segment(row, nu) = -du_max;
        u.segment(row, nu) = du_max;
        row += nu;
    }
}

template <typename Scalar, std::size_t NX>
inline void set_ellipsoidal_bounds(Eigen::VectorX<Scalar>& l, Eigen::VectorX<Scalar>& u, int& row, const ellipsoidal_set<Scalar, NX>& s)
{
    constexpr int nx = static_cast<int>(NX);
    Eigen::SelfAdjointEigenSolver<Matrix<Scalar, NX, NX>> eig(s.P);
    auto D = eig.eigenvalues();
    for(int i = 0; i < nx; ++i)
    {
        Scalar bound = std::sqrt(s.alpha / D(i));
        l(row + 2 * i) = -bound;
        u(row + 2 * i) = bound;
        l(row + 2 * i + 1) = -bound;
        u(row + 2 * i + 1) = bound;
    }
    row += 2 * nx;
}

template <typename Scalar, std::size_t NX>
inline void set_polytopic_bounds(Eigen::VectorX<Scalar>& l, Eigen::VectorX<Scalar>& u, int& row, const polytopic_set<Scalar, NX>& s)
{
    constexpr Scalar large_bound = Scalar{1e20};
    int n_faces = static_cast<int>(s.H.rows());
    l.segment(row, n_faces).setConstant(-large_bound);
    u.segment(row, n_faces) = s.h;
    row += n_faces;
}

template <typename Scalar, std::size_t NX>
inline void set_terminal_bounds(Eigen::VectorX<Scalar>& l, Eigen::VectorX<Scalar>& u, int& row, const terminal_set<Scalar, NX>& tset)
{
    std::visit([&](const auto& s) {
        using T = std::decay_t<decltype(s)>;
        if constexpr(std::is_same_v<T, ellipsoidal_set<Scalar, NX>>)
            set_ellipsoidal_bounds<Scalar, NX>(l, u, row, s);
        else
            set_polytopic_bounds<Scalar, NX>(l, u, row, s);
    }, tset);
}

template <typename Scalar, std::size_t NX, std::size_t NU>
[[nodiscard]] auto build_bounds_vectors(int N,
                                        const Vector<Scalar, NX>& x0,
                                        const std::optional<Vector<Scalar, NX>>& x_min,
                                        const std::optional<Vector<Scalar, NX>>& x_max,
                                        const std::optional<Vector<Scalar, NU>>& u_min,
                                        const std::optional<Vector<Scalar, NU>>& u_max,
                                        const std::optional<Vector<Scalar, NU>>& du_max,
                                        bool /*has_soft_constraints*/,
                                        const std::optional<terminal_set<Scalar, NX>>& terminal_constraint = {}) -> std::pair<Eigen::VectorX<Scalar>, Eigen::VectorX<Scalar>>
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);

    bool has_state_bounds = x_min.has_value() || x_max.has_value();
    bool has_input_bounds = u_min.has_value() || u_max.has_value();
    bool has_rate_bounds = du_max.has_value();

    int n_con = (N + 1) * nx + (has_state_bounds ? N * nx : 0) + (has_input_bounds ? N * nu : 0) + (has_rate_bounds ? N * nu : 0) + terminal_constraint_rows<Scalar, NX>(terminal_constraint);

    Eigen::VectorX<Scalar> l(n_con);
    Eigen::VectorX<Scalar> u(n_con);

    int row = 0;
    set_dynamics_bounds<Scalar, NX>(l, u, row, N, x0);
    if(has_state_bounds)
        set_state_bounds<Scalar, NX>(l, u, row, N, x_min, x_max);
    if(has_input_bounds)
        set_input_bounds<Scalar, NU>(l, u, row, N, u_min, u_max);
    if(has_rate_bounds)
        set_rate_bounds<Scalar, NU>(l, u, row, N, du_max.value());
    if(terminal_constraint.has_value())
        set_terminal_bounds<Scalar, NX>(l, u, row, terminal_constraint.value());

    return {l, u};
}

// --- Cost vector building ---

template <typename Scalar, std::size_t NX, std::size_t NU>
[[nodiscard]] auto build_cost_vector(int /*N*/, int n_dec, const Matrix<Scalar, NX, NX>& /*Q*/, const Matrix<Scalar, NX, NX>& /*Qf*/) -> Eigen::VectorX<Scalar>
{
    return Eigen::VectorX<Scalar>::Zero(n_dec);
}

template <typename Scalar, std::size_t NX, std::size_t NU>
[[nodiscard]] auto build_cost_vector(int N, int n_dec, const Matrix<Scalar, NX, NX>& Q, const Matrix<Scalar, NX, NX>& Qf, const Vector<Scalar, NX>& x_ref) -> Eigen::VectorX<Scalar>
{
    constexpr int nx = static_cast<int>(NX);
    Eigen::VectorX<Scalar> q = Eigen::VectorX<Scalar>::Zero(n_dec);
    Vector<Scalar, NX> Qxr = Q * x_ref;
    for(int k = 0; k < N; ++k)
        q.segment(k * nx, nx) = -Qxr;
    q.segment(N * nx, nx) = -(Qf * x_ref);
    return q;
}

template <typename Scalar, std::size_t NX, std::size_t NU>
[[nodiscard]] auto build_cost_vector(int N, int n_dec, const Matrix<Scalar, NX, NX>& Q, const Matrix<Scalar, NX, NX>& Qf, std::span<const Vector<Scalar, NX>> x_ref) -> Eigen::VectorX<Scalar>
{
    constexpr int nx = static_cast<int>(NX);
    Eigen::VectorX<Scalar> q = Eigen::VectorX<Scalar>::Zero(n_dec);
    for(int k = 0; k < N; ++k)
        q.segment(k * nx, nx) = -(Q * x_ref[static_cast<std::size_t>(k)]);
    q.segment(N * nx, nx) = -(Qf * x_ref[static_cast<std::size_t>(N)]);
    return q;
}

}

#endif
