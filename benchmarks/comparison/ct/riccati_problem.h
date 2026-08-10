#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_CT_RICCATI_PROBLEM_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_CT_RICCATI_PROBLEM_H

#include <Eigen/Dense>

#include <cstddef>
#include <algorithm>

namespace ctrlpp::bench
{

template <std::size_t NX, std::size_t NU>
struct riccati_plant
{
    Eigen::Matrix<double, int(NX), int(NX)> A;
    Eigen::Matrix<double, int(NX), int(NU)> B;
    Eigen::Matrix<double, int(NX), int(NX)> Q;
    Eigen::Matrix<double, int(NU), int(NU)> R;
};

// The spectrum is Re(lambda) < 0 on every rung, so the continuous Riccati
// equation is well-defined across the whole size sweep.
template <std::size_t NX, std::size_t NU>
riccati_plant<NX, NU> build_damped_chain()
{
    Eigen::Matrix<double, int(NX), int(NX)> A = Eigen::Matrix<double, int(NX), int(NX)>::Zero();
    for(std::size_t i = 0; i < NX; ++i)
        A(int(i), int(i)) = -0.5;
    for(std::size_t i = 0; i + 1 < NX; ++i)
        A(int(i), int(i + 1)) = 1.0;

    Eigen::Matrix<double, int(NX), int(NU)> B = Eigen::Matrix<double, int(NX), int(NU)>::Zero();
    const std::size_t group = NX / NU;
    for(std::size_t j = 0; j < NU; ++j)
    {
        const std::size_t last_row = std::min((j + 1) * group, NX) - 1;
        B(int(last_row), int(j)) = 1.0;
    }

    Eigen::Matrix<double, int(NX), int(NX)> Q = Eigen::Matrix<double, int(NX), int(NX)>::Identity();
    Eigen::Matrix<double, int(NU), int(NU)> R = 0.1 * Eigen::Matrix<double, int(NU), int(NU)>::Identity();
    return riccati_plant<NX, NU>{A, B, Q, R};
}

// A chain of integrators sampled at dt, so A is a transition matrix and the
// discrete Riccati equation is the one this plant poses.
template <std::size_t NX, std::size_t NU>
riccati_plant<NX, NU> build_integrator_chain(double dt)
{
    Eigen::Matrix<double, int(NX), int(NX)> A = Eigen::Matrix<double, int(NX), int(NX)>::Identity();
    for(std::size_t i = 0; i + 1 < NX; ++i)
        A(int(i), int(i + 1)) = dt;

    Eigen::Matrix<double, int(NX), int(NU)> B = Eigen::Matrix<double, int(NX), int(NU)>::Zero();
    const std::size_t group = NX / NU;
    for(std::size_t j = 0; j < NU; ++j)
    {
        const std::size_t last_row = std::min((j + 1) * group, NX) - 1;
        B(int(last_row), int(j)) = dt;
    }

    Eigen::Matrix<double, int(NX), int(NX)> Q = Eigen::Matrix<double, int(NX), int(NX)>::Identity();
    Eigen::Matrix<double, int(NU), int(NU)> R = 0.1 * Eigen::Matrix<double, int(NU), int(NU)>::Identity();
    return riccati_plant<NX, NU>{A, B, Q, R};
}

// The residual norm over the sum of the norms of the terms that cancel to form
// it, so the figure is dimensionless and comparable across the size sweep.
template <std::size_t NX, std::size_t NU>
double riccati_relative_residual(const riccati_plant<NX, NU>& plant, const Eigen::Matrix<double, int(NX), int(NX)>& P)
{
    const Eigen::Matrix<double, int(NX), int(NX)> cross = plant.A.transpose() * P + P * plant.A;
    const Eigen::Matrix<double, int(NX), int(NX)> quad =
        P * plant.B * plant.R.inverse() * plant.B.transpose() * P;
    return (cross - quad + plant.Q).norm() / (cross.norm() + quad.norm() + plant.Q.norm());
}

template <std::size_t NX, std::size_t NU>
double dare_relative_residual(const riccati_plant<NX, NU>& plant, const Eigen::Matrix<double, int(NX), int(NX)>& P)
{
    const Eigen::Matrix<double, int(NX), int(NX)> transported = plant.A.transpose() * P * plant.A;
    const Eigen::Matrix<double, int(NU), int(NU)> inner = plant.R + plant.B.transpose() * P * plant.B;
    const Eigen::Matrix<double, int(NU), int(NX)> cross = plant.B.transpose() * P * plant.A;
    const Eigen::Matrix<double, int(NX), int(NX)> quad = cross.transpose() * inner.inverse() * cross;
    return (transported - P - quad + plant.Q).norm()
           / (transported.norm() + P.norm() + quad.norm() + plant.Q.norm());
}

template <std::size_t NX, std::size_t NU>
double closed_loop_abscissa(const riccati_plant<NX, NU>& plant, const Eigen::Matrix<double, int(NU), int(NX)>& K)
{
    const Eigen::Matrix<double, int(NX), int(NX)> closed_loop = plant.A - plant.B * K;
    const Eigen::EigenSolver<Eigen::Matrix<double, int(NX), int(NX)>> spectrum(closed_loop, false);
    return spectrum.eigenvalues().real().maxCoeff();
}

template <std::size_t NX, std::size_t NU>
double closed_loop_spectral_radius(const riccati_plant<NX, NU>& plant, const Eigen::Matrix<double, int(NU), int(NX)>& K)
{
    const Eigen::Matrix<double, int(NX), int(NX)> closed_loop = plant.A - plant.B * K;
    const Eigen::EigenSolver<Eigen::Matrix<double, int(NX), int(NX)>> spectrum(closed_loop, false);
    return spectrum.eigenvalues().cwiseAbs().maxCoeff();
}

// Matrix form of X -> Ac' X + X Ac, laid out for Eigen's column-major vec:
// vec(Ac' X) = (I kron Ac') vec(X) and vec(X Ac) = (Ac' kron I) vec(X).
template <std::size_t NX>
Eigen::MatrixXd lyapunov_operator(const Eigen::Matrix<double, int(NX), int(NX)>& Ac)
{
    constexpr int n = int(NX);
    Eigen::MatrixXd op = Eigen::MatrixXd::Zero(n * n, n * n);
    for(int i = 0; i < n; ++i)
        for(int k = 0; k < n; ++k)
            for(int m = 0; m < n; ++m)
            {
                op(i * n + k, i * n + m) += Ac(m, k);
                op(i * n + k, m * n + k) += Ac(m, i);
            }
    return op;
}

// The cost matrix a gain actually realizes, recovered from that gain rather
// than from a second Riccati solve. For the optimal gain this is the Riccati
// solution the gain was formed from; for any other gain it is the closed-loop
// cost that gain achieves, whose Riccati residual is then the gain's own
// distance from optimality.
template <std::size_t NX, std::size_t NU>
Eigen::Matrix<double, int(NX), int(NX)> realized_cost_matrix(const riccati_plant<NX, NU>& plant,
                                                             const Eigen::Matrix<double, int(NU), int(NX)>& K)
{
    constexpr int n = int(NX);
    const Eigen::Matrix<double, n, n> closed_loop = plant.A - plant.B * K;
    const Eigen::Matrix<double, n, n> cost = plant.Q + K.transpose() * plant.R * K;
    const Eigen::VectorXd rhs = -Eigen::Map<const Eigen::VectorXd>(cost.data(), n * n);
    const Eigen::VectorXd solution = lyapunov_operator<NX>(closed_loop).colPivHouseholderQr().solve(rhs);
    return Eigen::Map<const Eigen::Matrix<double, n, n>>(solution.data());
}

// Matrix form of X -> Ac' X Ac, laid out for Eigen's column-major vec:
// vec(Ac' X Ac) = (Ac' kron Ac') vec(X).
template <std::size_t NX>
Eigen::MatrixXd stein_operator(const Eigen::Matrix<double, int(NX), int(NX)>& Ac)
{
    constexpr int n = int(NX);
    Eigen::MatrixXd op(n * n, n * n);
    for(int i = 0; i < n; ++i)
        for(int k = 0; k < n; ++k)
            for(int j = 0; j < n; ++j)
                for(int m = 0; m < n; ++m)
                    op(i * n + k, j * n + m) = Ac(j, i) * Ac(m, k);
    return op;
}

// The discrete counterpart of realized_cost_matrix: the cost a gain achieves on
// a sampled plant, recovered from that gain through the closed-loop Stein
// equation rather than from a second Riccati solve.
template <std::size_t NX, std::size_t NU>
Eigen::Matrix<double, int(NX), int(NX)> realized_cost_matrix_discrete(
    const riccati_plant<NX, NU>& plant, const Eigen::Matrix<double, int(NU), int(NX)>& K)
{
    constexpr int n = int(NX);
    const Eigen::Matrix<double, n, n> closed_loop = plant.A - plant.B * K;
    const Eigen::Matrix<double, n, n> cost = plant.Q + K.transpose() * plant.R * K;
    const Eigen::VectorXd rhs = -Eigen::Map<const Eigen::VectorXd>(cost.data(), n * n);
    const Eigen::MatrixXd op = stein_operator<NX>(closed_loop) - Eigen::MatrixXd::Identity(n * n, n * n);
    const Eigen::VectorXd solution = op.colPivHouseholderQr().solve(rhs);
    return Eigen::Map<const Eigen::Matrix<double, n, n>>(solution.data());
}

template <std::size_t NX, std::size_t NU>
double discrete_gain_optimality_residual(const riccati_plant<NX, NU>& plant,
                                         const Eigen::Matrix<double, int(NU), int(NX)>& K)
{
    return dare_relative_residual<NX, NU>(plant, realized_cost_matrix_discrete<NX, NU>(plant, K));
}

// How far a gain is from being optimal for the cost it realizes: zero exactly
// when K = R^-1 B' P for that gain's own P, and quadratic in the gain error
// elsewhere, so it is an upper bound rather than a sensitive discriminator.
template <std::size_t NX, std::size_t NU>
double gain_optimality_residual(const riccati_plant<NX, NU>& plant, const Eigen::Matrix<double, int(NU), int(NX)>& K)
{
    return riccati_relative_residual<NX, NU>(plant, realized_cost_matrix<NX, NU>(plant, K));
}

}

#endif
