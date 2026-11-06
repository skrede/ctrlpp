#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_PREDICTIVE_BOUND_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_PREDICTIVE_BOUND_H

#include "golden_reference.h"
#include "derived_tolerance.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <algorithm>

namespace ctrlpp
{

// A step of the backward Riccati recursion forms P*B, B^T*P*B + R, P*A and B^T*P*A, solves for K, then forms
// A - B*K and Q + A^T*P*(A - B*K). Products are charged 2*p*q*r, sums one per entry and the solve
// 4*nu^3 + 2*nu^2*nx, as the Riccati stages are; like gain_departure_bound, the count treats the recursion
// as forward-stable to it and carries no condition number.
constexpr double finite_horizon_roundings(std::size_t nx, std::size_t nu, std::size_t horizon)
{
    const double n = static_cast<double>(nx);
    const double p = static_cast<double>(nu);
    const double step = 6.0 * n * n * n + 6.0 * n * n * p + 4.0 * p * p * n + 4.0 * p * p * p + p * p + 2.0 * n * n;
    return static_cast<double>(horizon) * step + 2.0 * n * n + n * p + p * p;
}

constexpr double predictive_gain_departure()
{
    const double gain_norm = (golden::host_predictive_k0 < 0.0 ? -golden::host_predictive_k0 : golden::host_predictive_k0) +
            (golden::host_predictive_k1 < 0.0 ? -golden::host_predictive_k1 : golden::host_predictive_k1);
    return counted_departure_bound<double>(finite_horizon_roundings(golden::predictive_nx, golden::predictive_nu, golden::predictive_nh), gain_norm);
}

// An upper bound on the l1 norm of the inverse KKT matrix's first-input row. The host solve's relative
// error is at most the condition number times its operation count, a 4*m^3 LU at m = NV + MaxM.
constexpr double predictive_input_row()
{
    const double m = static_cast<double>(golden::predictive_nv + golden::predictive_max_m);
    return golden::host_predictive_input_row * (1.0 + golden::host_predictive_kkt_condition * counted_departure_bound<double>(4.0 * m * m * m, 1.0));
}

// A lower bound on the constraint Jacobian's smallest singular value. Singular values move by at most the
// perturbation's 2-norm, and the host SVD is charged 4*m^3 roundings at the Jacobian's own norm.
constexpr double predictive_jacobian_sigma()
{
    const double m = static_cast<double>(golden::predictive_nv);
    return golden::host_predictive_jacobian_sigma - counted_departure_bound<double>(4.0 * m * m * m, golden::host_predictive_jacobian_norm);
}

// The largest error of a finite-differenced dynamics entry of the constraint Jacobian, for a decision
// vector no larger than z_inf. Each central difference evaluates the continuity row z_next - f(x, u) twice,
// three roundings each at S = (2 + dt)*(z_inf + h), subtracts (charged at 2*S) and divides by 2*h. The step
// is cbrt(eps)*max(1, |z_j|) made exactly representable, so within a factor of two of that: the lower end
// sets the divisor and the upper end S. The evaluation point z_j - h rounds once, moving the exact quotient
// of a unit-slope row by eps*(z_inf + h)/(2*h), and the quotient itself rounds once more.
inline double fd_entry_bound(double z_inf)
{
    const double eps = std::numeric_limits<double>::epsilon();
    const double root = std::cbrt(eps);
    const double h_low = root / 2.0;
    const double h_high = 2.0 * root * std::max(1.0, z_inf);
    const double scale = (2.0 + golden::predictive_dt) * (z_inf + h_high);
    return 8.0 * eps * scale / (2.0 * h_low) + eps * (z_inf + h_high) / (2.0 * h_low) + 2.0 * eps;
}

// From the solver's stationarity leg: any lambda with |H*z - J~^T*lambda| <= s is (J~*J~^T)^-1*J~*(H*z - s),
// whose 2-norm is at most sqrt(nv)*(|H*z| + s) over the smallest singular value of J~, the
// finite-differenced Jacobian.
inline double multiplier_bound(double z_inf, double kkt_tol, double stationarity_charge, double charge_slope)
{
    const double nv = static_cast<double>(golden::predictive_nv);
    const double fd = fd_entry_bound(z_inf);
    const double spread = std::sqrt(static_cast<double>(golden::predictive_nx * (golden::predictive_nx + golden::predictive_nu))) * fd;
    const double sigma = predictive_jacobian_sigma() - spread;
    const double h_norm = std::max(golden::predictive_state_weight, golden::predictive_input_weight);
    const double numer = std::sqrt(nv) * (h_norm * z_inf + kkt_tol + stationarity_charge);
    return numer / (sigma - std::sqrt(nv) * charge_slope);
}

// The solver measures stationarity with the finite-differenced Jacobian, so the exact residual adds
// |(J~ - J)^T*lambda|, at most nx entry errors per column. Forming the gradient, J~^T*lambda and their
// difference rounds 2*nx + 4 times; each feasibility entry rounds three times at (2 + dt)*z_inf.
inline double kkt_residual_bound(double z_inf, double kkt_tol)
{
    const double eps = std::numeric_limits<double>::epsilon();
    const double roundings = (2.0 * static_cast<double>(golden::predictive_nx) + 4.0) * eps;
    const double h_norm = std::max(golden::predictive_state_weight, golden::predictive_input_weight);
    const double fd = fd_entry_bound(z_inf);
    const double slope = roundings * static_cast<double>(golden::predictive_nx + 1) * (1.0 + fd);
    const double lambda = multiplier_bound(z_inf, kkt_tol, roundings * h_norm * z_inf, slope);
    const double stationary = kkt_tol + roundings * h_norm * z_inf + slope * lambda + static_cast<double>(golden::predictive_nx) * fd * lambda;
    const double feasibility = kkt_tol + counted_departure_bound<double>(3.0, (2.0 + golden::predictive_dt) * z_inf);
    return std::max(stationary, feasibility);
}

// The posed problem is an equality-constrained QP, so for any (z, lambda) the exact KKT matrix M maps
// (z - z*, lambda - lambda*) to the exact residual, and |u0 - u0*| is at most the l1 norm of M^-1's
// first-input row times the residual's infinity norm.
// @cite nocedal2006 -- Nocedal & Wright, "Numerical Optimization", 2nd ed., 2006, Sec. 16.1 (the KKT system of an equality-constrained QP)
inline double input_certificate(double z_inf, double kkt_tol)
{
    return predictive_input_row() * kkt_residual_bound(z_inf, kkt_tol);
}

}

#endif
