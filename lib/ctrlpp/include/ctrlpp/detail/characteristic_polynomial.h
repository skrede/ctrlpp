#ifndef HPP_GUARD_CTRLPP_DETAIL_CHARACTERISTIC_POLYNOMIAL_H
#define HPP_GUARD_CTRLPP_DETAIL_CHARACTERISTIC_POLYNOMIAL_H

/// @brief Characteristic polynomial of a real square matrix from the eigenvalues of its balanced copy.
///
/// The balance is an exact power-of-two similarity, so the eigenvalues are the matrix's own; the QR
/// algorithm computes them backward stably, and the monic product of their real linear and quadratic
/// factors forms the coefficients, highest degree first.
///
/// @cite golub2013 -- Golub & Van Loan, "Matrix Computations", 4th ed., 2013, Sec. 7.5 (QR algorithm)

#include "ctrlpp/expected.h"

#include "ctrlpp/detail/diagonal_balance.h"

#include "ctrlpp/model/conversion_types.h"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <array>
#include <complex>
#include <cstddef>

namespace ctrlpp::detail
{

/// Real eigenvalues contribute (s - lambda); a conjugate pair a +- b i, adjacent in Eigen's output,
/// contributes s^2 - 2 a s + (a^2 + b^2) once. Coefficients above the running degree are zero, so each
/// factor sweeps the whole array without changing them.
template <typename Scalar, int N>
auto monic_from_eigenvalues(const Eigen::Matrix<std::complex<Scalar>, N, 1>& lambda) -> std::array<Scalar, std::size_t(N) + 1>
{
    constexpr auto n = static_cast<std::size_t>(N);
    std::array<Scalar, n + 1> p{};
    p[0] = Scalar{1};
    for(int i = 0; i < N; ++i)
    {
        const Scalar a = lambda(i).real();
        const Scalar b = lambda(i).imag();
        if(b == Scalar{0})
        {
            for(std::size_t k = n; k >= 1; --k)
                p[k] -= a * p[k - 1];
            continue;
        }
        const Scalar c1 = -Scalar{2} * a;
        const Scalar c0 = a * a + b * b;
        for(std::size_t k = n; k >= 2; --k)
            p[k] += c1 * p[k - 1] + c0 * p[k - 2];
        p[1] += c1 * p[0];
        ++i;
    }
    return p;
}

template <typename Scalar, int N>
auto characteristic_polynomial(Eigen::Matrix<Scalar, N, N> M)
    -> ctrlpp::expected<std::array<Scalar, std::size_t(N) + 1>, conversion_error>
{
    if(!M.allFinite())
        return ctrlpp::unexpected(conversion_error::unrepresentable_result);
    Eigen::Matrix<Scalar, N, 1> scaling;
    balance_diagonal<Scalar, N>(M, scaling);
    const Eigen::EigenSolver<Eigen::Matrix<Scalar, N, N>> es(M, false);
    if(es.info() != Eigen::Success)
        return ctrlpp::unexpected(conversion_error::eigenvalues_not_converged);
    return monic_from_eigenvalues<Scalar, N>(es.eigenvalues());
}

}

#endif
