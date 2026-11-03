#ifndef HPP_GUARD_CTRLPP_DETAIL_QUASI_TRIANGULAR_H
#define HPP_GUARD_CTRLPP_DETAIL_QUASI_TRIANGULAR_H

// A 2x2 block's eigenvalues are the half-trace plus and minus the square root of the discriminant, not its
// diagonal entries: those are neither root of a real pair, and the shared real part of a conjugate pair
// only once the block is standardized.
// @cite golub2013 -- Golub & Van Loan, "Matrix Computations", 4th ed., 2013, Sec. 7.4.1
// @cite bai_demmel_1993 -- Bai & Demmel, "On swapping diagonal blocks in real Schur form", 1993

#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <complex>
#include <algorithm>

namespace ctrlpp::detail
{

// For a 2x2 block `first` is the root with the larger real part, so a half-plane test reads `first` alone.
template <typename Scalar>
struct block_spectrum
{
    std::complex<Scalar> first;
    std::complex<Scalar> second;
    bool real_pair;
};

// The LAPACK significance test on the computed subdiagonal entry: unit roundoff times the larger of the
// whole factor's magnitude and the block's own diagonal magnitude. `factor_scale` is passed in so a walk
// over every block stays linear in the factor's size.
template <typename Scalar, int N>
Eigen::Index quasi_triangular_block_size(const Eigen::Matrix<Scalar, N, N>& T, Eigen::Index pos, Scalar factor_scale)
{
    using std::abs;

    if(pos + 1 >= N)
        return 1;

    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    const Scalar local_scale = std::max(abs(T(pos, pos)), abs(T(pos + 1, pos + 1)));
    return (abs(T(pos + 1, pos)) > eps * std::max(factor_scale, local_scale)) ? 2 : 1;
}

template <typename Scalar>
block_spectrum<Scalar> block_spectrum_from_discriminant(Scalar half_trace, Scalar discriminant)
{
    using std::sqrt;

    if(discriminant >= Scalar{0})
    {
        const Scalar root = sqrt(discriminant);
        return {std::complex<Scalar>(half_trace + root, Scalar{0}), std::complex<Scalar>(half_trace - root, Scalar{0}), true};
    }

    const Scalar imaginary_part = sqrt(-discriminant);
    return {std::complex<Scalar>(half_trace, imaginary_part), std::complex<Scalar>(half_trace, -imaginary_part), false};
}

template <typename Scalar, int N>
block_spectrum<Scalar> quasi_triangular_block_spectrum(const Eigen::Matrix<Scalar, N, N>& T, Eigen::Index pos, Eigen::Index block_size)
{
    if(block_size == 1)
    {
        const std::complex<Scalar> eigenvalue(T(pos, pos), Scalar{0});
        return {eigenvalue, eigenvalue, true};
    }

    const Scalar a = T(pos, pos);
    const Scalar b = T(pos, pos + 1);
    const Scalar c = T(pos + 1, pos);
    const Scalar d = T(pos + 1, pos + 1);
    const Scalar half_trace = (a + d) / Scalar{2};
    const Scalar determinant = a * d - b * c;
    const Scalar discriminant = half_trace * half_trace - determinant;
    return block_spectrum_from_discriminant(half_trace, discriminant);
}

// A backward-stable factorization returns the exact factor of a nearby matrix, so an eigenvalue read from
// an N-by-N factor carries a perturbation of order N times unit roundoff times the factor's magnitude.
template <typename Scalar, int N>
Scalar schur_eigenvalue_margin(const Eigen::Matrix<Scalar, N, N>& T)
{
    return Scalar{N} * std::numeric_limits<Scalar>::epsilon() * T.cwiseAbs().maxCoeff();
}

// Reading the half-trace instead of `first` would certify a 2x2 block whose two real roots have opposite
// signs, since its half-trace can sit left of the margin while its larger root is in the right half-plane.
template <typename Scalar, int N>
bool quasi_triangular_spectrum_strictly_left_of(const Eigen::Matrix<Scalar, N, N>& T, Scalar margin)
{
    const Scalar factor_scale = T.cwiseAbs().maxCoeff();

    Eigen::Index index = 0;
    while(index < N)
    {
        const Eigen::Index block_size = quasi_triangular_block_size<Scalar, N>(T, index, factor_scale);
        const block_spectrum<Scalar> block = quasi_triangular_block_spectrum<Scalar, N>(T, index, block_size);
        if(!(block.first.real() < -margin))
            return false;
        index += block_size;
    }
    return true;
}

}

#endif
