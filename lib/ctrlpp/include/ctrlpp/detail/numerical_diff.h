#ifndef HPP_GUARD_CTRLPP_DETAIL_NUMERICAL_DIFF_H
#define HPP_GUARD_CTRLPP_DETAIL_NUMERICAL_DIFF_H

#include "ctrlpp/types.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <span>
#include <vector>

namespace ctrlpp::detail
{

// ---------------------------------------------------------------------------
// Span-based finite-difference utilities (NLopt compatibility)
// ---------------------------------------------------------------------------

template <typename Scalar>
void finite_diff_gradient(const std::function<Scalar(std::span<const Scalar>)>& f, std::span<const Scalar> z, std::span<Scalar> grad)
{
    const auto eps = std::sqrt(std::numeric_limits<Scalar>::epsilon());
    const auto n = z.size();

    std::vector<Scalar> z_mut(z.begin(), z.end());

    for(std::size_t j = 0; j < n; ++j)
    {
        const Scalar h = eps * std::max(Scalar{1}, std::abs(z[j]));
        const Scalar orig = z_mut[j];

        z_mut[j] = orig + h;
        const Scalar f_plus = f(std::span<const Scalar>{z_mut.data(), n});

        z_mut[j] = orig - h;
        const Scalar f_minus = f(std::span<const Scalar>{z_mut.data(), n});

        grad[j] = (f_plus - f_minus) / (Scalar{2} * h);
        z_mut[j] = orig;
    }
}

template <typename Scalar>
void finite_diff_jacobian(const std::function<void(std::span<const Scalar>, std::span<Scalar>)>& c, int n_constraints, std::span<const Scalar> z, std::span<Scalar> jac)
{
    const auto eps = std::sqrt(std::numeric_limits<Scalar>::epsilon());
    const auto n = z.size();
    const auto m = static_cast<std::size_t>(n_constraints);

    std::vector<Scalar> z_mut(z.begin(), z.end());
    std::vector<Scalar> c_plus(m);
    std::vector<Scalar> c_minus(m);

    for(std::size_t j = 0; j < n; ++j)
    {
        const Scalar h = eps * std::max(Scalar{1}, std::abs(z[j]));
        const Scalar orig = z_mut[j];

        z_mut[j] = orig + h;
        c(std::span<const Scalar>{z_mut.data(), n}, std::span<Scalar>{c_plus.data(), m});

        z_mut[j] = orig - h;
        c(std::span<const Scalar>{z_mut.data(), n}, std::span<Scalar>{c_minus.data(), m});

        for(std::size_t i = 0; i < m; ++i)
        {
            jac[i * n + j] = (c_plus[i] - c_minus[i]) / (Scalar{2} * h);
        }

        z_mut[j] = orig;
    }
}

// ---------------------------------------------------------------------------
// Eigen-native numerical Jacobian utilities (zero heap allocation)
// ---------------------------------------------------------------------------

/// Compute dF/dx via central differences for f(x, u) -> x_next.
/// Returns a fixed-size NX x NX Jacobian matrix.
template <typename Scalar, std::size_t NX, std::size_t NU, typename F>
auto numerical_jacobian_x(const F& f, const Vector<Scalar, NX>& x, const Vector<Scalar, NU>& u, Scalar eps = std::sqrt(std::numeric_limits<Scalar>::epsilon())) -> Matrix<Scalar, NX, NX>
{
    Matrix<Scalar, NX, NX> jac;
    Vector<Scalar, NX> x_perturbed = x;

    for(std::size_t j = 0; j < NX; ++j)
    {
        const auto idx = static_cast<Eigen::Index>(j);
        const Scalar h = eps * std::max(Scalar{1}, std::abs(x[idx]));

        x_perturbed[idx] = x[idx] + h;
        const Vector<Scalar, NX> f_plus = f(x_perturbed, u);

        x_perturbed[idx] = x[idx] - h;
        const Vector<Scalar, NX> f_minus = f(x_perturbed, u);

        jac.col(idx) = (f_plus - f_minus) / (Scalar{2} * h);

        x_perturbed[idx] = x[idx];
    }

    return jac;
}

/// Compute dF/du via central differences for f(x, u) -> x_next.
/// Returns a fixed-size NX x NU Jacobian matrix.
template <typename Scalar, std::size_t NX, std::size_t NU, typename F>
auto numerical_jacobian_u(const F& f, const Vector<Scalar, NX>& x, const Vector<Scalar, NU>& u, Scalar eps = std::sqrt(std::numeric_limits<Scalar>::epsilon())) -> Matrix<Scalar, NX, NU>
{
    Matrix<Scalar, NX, NU> jac;
    Vector<Scalar, NU> u_perturbed = u;

    for(std::size_t j = 0; j < NU; ++j)
    {
        const auto idx = static_cast<Eigen::Index>(j);
        const Scalar h = eps * std::max(Scalar{1}, std::abs(u[idx]));

        u_perturbed[idx] = u[idx] + h;
        const Vector<Scalar, NX> f_plus = f(x, u_perturbed);

        u_perturbed[idx] = u[idx] - h;
        const Vector<Scalar, NX> f_minus = f(x, u_perturbed);

        jac.col(idx) = (f_plus - f_minus) / (Scalar{2} * h);

        u_perturbed[idx] = u[idx];
    }

    return jac;
}

/// Compute dH/dx via central differences for h(x) -> y.
/// Returns a fixed-size NY x NX Jacobian matrix.
template <typename Scalar, std::size_t NX, std::size_t NY, typename H>
auto numerical_jacobian_h(const H& h, const Vector<Scalar, NX>& x, Scalar eps = std::sqrt(std::numeric_limits<Scalar>::epsilon())) -> Matrix<Scalar, NY, NX>
{
    Matrix<Scalar, NY, NX> jac;
    Vector<Scalar, NX> x_perturbed = x;

    for(std::size_t j = 0; j < NX; ++j)
    {
        const auto idx = static_cast<Eigen::Index>(j);
        const Scalar step = eps * std::max(Scalar{1}, std::abs(x[idx]));

        x_perturbed[idx] = x[idx] + step;
        const Vector<Scalar, NY> h_plus = h(x_perturbed);

        x_perturbed[idx] = x[idx] - step;
        const Vector<Scalar, NY> h_minus = h(x_perturbed);

        jac.col(idx) = (h_plus - h_minus) / (Scalar{2} * step);

        x_perturbed[idx] = x[idx];
    }

    return jac;
}

/// Compute dg/dx via central differences for g(x, u) -> Vector<NC>.
/// Returns a fixed-size NC x NX Jacobian sub-block.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC, typename G>
auto numerical_jacobian_gx(const G& g, const Vector<Scalar, NX>& x, const Vector<Scalar, NU>& u, Scalar eps = std::sqrt(std::numeric_limits<Scalar>::epsilon())) -> Matrix<Scalar, NC, NX>
{
    Matrix<Scalar, NC, NX> jac;
    Vector<Scalar, NX> x_perturbed = x;

    for(std::size_t j = 0; j < NX; ++j)
    {
        const auto idx = static_cast<Eigen::Index>(j);
        const Scalar h = eps * std::max(Scalar{1}, std::abs(x[idx]));

        x_perturbed[idx] = x[idx] + h;
        const Vector<Scalar, NC> g_plus = g(x_perturbed, u);

        x_perturbed[idx] = x[idx] - h;
        const Vector<Scalar, NC> g_minus = g(x_perturbed, u);

        jac.col(idx) = (g_plus - g_minus) / (Scalar{2} * h);
        x_perturbed[idx] = x[idx];
    }

    return jac;
}

/// Compute dg/du via central differences for g(x, u) -> Vector<NC>.
/// Returns a fixed-size NC x NU Jacobian sub-block.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC, typename G>
auto numerical_jacobian_gu(const G& g, const Vector<Scalar, NX>& x, const Vector<Scalar, NU>& u, Scalar eps = std::sqrt(std::numeric_limits<Scalar>::epsilon())) -> Matrix<Scalar, NC, NU>
{
    Matrix<Scalar, NC, NU> jac;
    Vector<Scalar, NU> u_perturbed = u;

    for(std::size_t j = 0; j < NU; ++j)
    {
        const auto idx = static_cast<Eigen::Index>(j);
        const Scalar h = eps * std::max(Scalar{1}, std::abs(u[idx]));

        u_perturbed[idx] = u[idx] + h;
        const Vector<Scalar, NC> g_plus = g(x, u_perturbed);

        u_perturbed[idx] = u[idx] - h;
        const Vector<Scalar, NC> g_minus = g(x, u_perturbed);

        jac.col(idx) = (g_plus - g_minus) / (Scalar{2} * h);
        u_perturbed[idx] = u[idx];
    }

    return jac;
}

/// Compute dg/d[x;u] via central differences for g(x, u) -> Vector<NC>.
/// Returns a fixed-size NC x (NX+NU) combined Jacobian matrix.
/// Columns 0..NX-1 correspond to dg/dx, columns NX..NX+NU-1 to dg/du.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC, typename G>
auto numerical_jacobian_g(const G& g, const Vector<Scalar, NX>& x, const Vector<Scalar, NU>& u, Scalar eps = std::sqrt(std::numeric_limits<Scalar>::epsilon())) -> Matrix<Scalar, NC, NX + NU>
{
    Matrix<Scalar, NC, NX + NU> jac;
    jac.template leftCols<NX>() = numerical_jacobian_gx<Scalar, NX, NU, NC>(g, x, u, eps);
    jac.template rightCols<NU>() = numerical_jacobian_gu<Scalar, NX, NU, NC>(g, x, u, eps);
    return jac;
}

/// Compute dh/dx via central differences for terminal constraint h(x) -> Vector<NTC>.
/// Returns a fixed-size NTC x NX Jacobian matrix.
template <typename Scalar, std::size_t NX, std::size_t NTC, typename H>
auto numerical_jacobian_tc(const H& h, const Vector<Scalar, NX>& x, Scalar eps = std::sqrt(std::numeric_limits<Scalar>::epsilon())) -> Matrix<Scalar, NTC, NX>
{
    return numerical_jacobian_h<Scalar, NX, NTC>(h, x, eps);
}

} // namespace ctrlpp::detail

#endif
