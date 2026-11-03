#ifndef HPP_GUARD_CTRLPP_DETAIL_NUMERICAL_JACOBIANS_H
#define HPP_GUARD_CTRLPP_DETAIL_NUMERICAL_JACOBIANS_H

// Central-difference Jacobians of fixed-size Eigen callables, allocation-free, with the step
// h = cbrt(eps) * max(1, |x_j|) rounded to a representable difference.
// @cite nocedal2006  -- Nocedal & Wright, "Numerical Optimization", 2nd ed., 2006, Sec. 8.1 (central-scheme FD step h = eps^(1/3)*scale; representable-h rounding p. 196)

#include "ctrlpp/types.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, typename F>
Matrix<Scalar, NX, NX> numerical_jacobian_x(const F& f, const Vector<Scalar, NX>& x, const Vector<Scalar, NU>& u, Scalar step_scale = std::cbrt(std::numeric_limits<Scalar>::epsilon()))
{
    Matrix<Scalar, NX, NX> jac;
    Vector<Scalar, NX> x_perturbed = x;

    for(std::size_t j = 0; j < NX; ++j)
    {
        const auto idx = static_cast<Eigen::Index>(j);
        const Scalar h_raw = step_scale * std::max(Scalar{1}, std::abs(x[idx]));
        const Scalar temp = x[idx] + h_raw;
        const Scalar h = temp - x[idx];

        x_perturbed[idx] = x[idx] + h;
        const Vector<Scalar, NX> f_plus = f(x_perturbed, u);

        x_perturbed[idx] = x[idx] - h;
        const Vector<Scalar, NX> f_minus = f(x_perturbed, u);

        jac.col(idx) = (f_plus - f_minus) / (Scalar{2} * h);

        x_perturbed[idx] = x[idx];
    }

    return jac;
}

template <typename Scalar, std::size_t NX, std::size_t NU, typename F>
Matrix<Scalar, NX, NU> numerical_jacobian_u(const F& f, const Vector<Scalar, NX>& x, const Vector<Scalar, NU>& u, Scalar step_scale = std::cbrt(std::numeric_limits<Scalar>::epsilon()))
{
    Matrix<Scalar, NX, NU> jac;
    Vector<Scalar, NU> u_perturbed = u;

    for(std::size_t j = 0; j < NU; ++j)
    {
        const auto idx = static_cast<Eigen::Index>(j);
        const Scalar h_raw = step_scale * std::max(Scalar{1}, std::abs(u[idx]));
        const Scalar temp = u[idx] + h_raw;
        const Scalar h = temp - u[idx];

        u_perturbed[idx] = u[idx] + h;
        const Vector<Scalar, NX> f_plus = f(x, u_perturbed);

        u_perturbed[idx] = u[idx] - h;
        const Vector<Scalar, NX> f_minus = f(x, u_perturbed);

        jac.col(idx) = (f_plus - f_minus) / (Scalar{2} * h);

        u_perturbed[idx] = u[idx];
    }

    return jac;
}

template <typename Scalar, std::size_t NX, std::size_t NY, typename H>
Matrix<Scalar, NY, NX> numerical_jacobian_h(const H& h, const Vector<Scalar, NX>& x, Scalar step_scale = std::cbrt(std::numeric_limits<Scalar>::epsilon()))
{
    Matrix<Scalar, NY, NX> jac;
    Vector<Scalar, NX> x_perturbed = x;

    for(std::size_t j = 0; j < NX; ++j)
    {
        const auto idx = static_cast<Eigen::Index>(j);
        const Scalar step_raw = step_scale * std::max(Scalar{1}, std::abs(x[idx]));
        const Scalar temp = x[idx] + step_raw;
        const Scalar step = temp - x[idx];

        x_perturbed[idx] = x[idx] + step;
        const Vector<Scalar, NY> h_plus = h(x_perturbed);

        x_perturbed[idx] = x[idx] - step;
        const Vector<Scalar, NY> h_minus = h(x_perturbed);

        jac.col(idx) = (h_plus - h_minus) / (Scalar{2} * step);

        x_perturbed[idx] = x[idx];
    }

    return jac;
}

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC, typename G>
Matrix<Scalar, NC, NX> numerical_jacobian_gx(const G& g, const Vector<Scalar, NX>& x, const Vector<Scalar, NU>& u, Scalar step_scale = std::cbrt(std::numeric_limits<Scalar>::epsilon()))
{
    Matrix<Scalar, NC, NX> jac;
    Vector<Scalar, NX> x_perturbed = x;

    for(std::size_t j = 0; j < NX; ++j)
    {
        const auto idx = static_cast<Eigen::Index>(j);
        const Scalar h_raw = step_scale * std::max(Scalar{1}, std::abs(x[idx]));
        const Scalar temp = x[idx] + h_raw;
        const Scalar h = temp - x[idx];

        x_perturbed[idx] = x[idx] + h;
        const Vector<Scalar, NC> g_plus = g(x_perturbed, u);

        x_perturbed[idx] = x[idx] - h;
        const Vector<Scalar, NC> g_minus = g(x_perturbed, u);

        jac.col(idx) = (g_plus - g_minus) / (Scalar{2} * h);
        x_perturbed[idx] = x[idx];
    }

    return jac;
}

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC, typename G>
Matrix<Scalar, NC, NU> numerical_jacobian_gu(const G& g, const Vector<Scalar, NX>& x, const Vector<Scalar, NU>& u, Scalar step_scale = std::cbrt(std::numeric_limits<Scalar>::epsilon()))
{
    Matrix<Scalar, NC, NU> jac;
    Vector<Scalar, NU> u_perturbed = u;

    for(std::size_t j = 0; j < NU; ++j)
    {
        const auto idx = static_cast<Eigen::Index>(j);
        const Scalar h_raw = step_scale * std::max(Scalar{1}, std::abs(u[idx]));
        const Scalar temp = u[idx] + h_raw;
        const Scalar h = temp - u[idx];

        u_perturbed[idx] = u[idx] + h;
        const Vector<Scalar, NC> g_plus = g(x, u_perturbed);

        u_perturbed[idx] = u[idx] - h;
        const Vector<Scalar, NC> g_minus = g(x, u_perturbed);

        jac.col(idx) = (g_plus - g_minus) / (Scalar{2} * h);
        u_perturbed[idx] = u[idx];
    }

    return jac;
}

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC, typename G>
Matrix<Scalar, NC, NX + NU> numerical_jacobian_g(
        const G& g, const Vector<Scalar, NX>& x, const Vector<Scalar, NU>& u, Scalar step_scale = std::cbrt(std::numeric_limits<Scalar>::epsilon()))
{
    Matrix<Scalar, NC, NX + NU> jac;
    jac.template leftCols<NX>() = numerical_jacobian_gx<Scalar, NX, NU, NC>(g, x, u, step_scale);
    jac.template rightCols<NU>() = numerical_jacobian_gu<Scalar, NX, NU, NC>(g, x, u, step_scale);
    return jac;
}

template <typename Scalar, std::size_t NX, std::size_t NTC, typename H>
Matrix<Scalar, NTC, NX> numerical_jacobian_tc(const H& h, const Vector<Scalar, NX>& x, Scalar step_scale = std::cbrt(std::numeric_limits<Scalar>::epsilon()))
{
    return numerical_jacobian_h<Scalar, NX, NTC>(h, x, step_scale);
}

}

#endif
