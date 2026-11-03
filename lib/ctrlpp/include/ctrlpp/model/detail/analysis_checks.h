#ifndef HPP_GUARD_CTRLPP_MODEL_DETAIL_ANALYSIS_CHECKS_H
#define HPP_GUARD_CTRLPP_MODEL_DETAIL_ANALYSIS_CHECKS_H

#include "ctrlpp/types.h"

#include <Eigen/Eigenvalues>

#include <array>
#include <cmath>
#include <limits>
#include <complex>
#include <cstddef>

namespace ctrlpp::detail
{

// A non-finite matrix and a solver result other than Eigen::Success both fill the whole array with a quiet
// NaN, so an unspecified or partially populated eigenvalue vector is never mistaken for a spectrum.
template <typename Scalar, std::size_t NX>
std::array<std::complex<Scalar>, NX> spectrum(const Matrix<Scalar, NX, NX>& A)
{
    constexpr int n = static_cast<int>(NX);
    constexpr Scalar nan = std::numeric_limits<Scalar>::quiet_NaN();

    std::array<std::complex<Scalar>, NX> result;
    result.fill(std::complex<Scalar>{nan, nan});

    if(!A.allFinite())
        return result;

    Eigen::EigenSolver<Eigen::Matrix<Scalar, n, n>> solver(A, false);
    if(solver.info() != Eigen::Success)
        return result;

    auto evals = solver.eigenvalues();
    for(Eigen::Index i = 0; i < n; ++i)
        result[static_cast<std::size_t>(i)] = evals(i);
    return result;
}

template <typename Scalar>
bool is_determinate(const std::complex<Scalar>& pole)
{
    return std::isfinite(pole.real()) && std::isfinite(pole.imag());
}

}

#endif
