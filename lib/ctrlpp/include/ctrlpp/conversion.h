#ifndef HPP_GUARD_CTRLPP_CONVERSION_H
#define HPP_GUARD_CTRLPP_CONVERSION_H

#include "ctrlpp/state_space.h"
#include "ctrlpp/transfer_function.h"

#include <array>
#include <cstddef>

namespace ctrlpp {

// tf2ss: Transfer function to controllable canonical form state-space.
// H(s) = num(s) / den(s), coefficients highest-degree-first (MATLAB convention).
// Requires NumDeg <= DenDeg (proper transfer function).
// Returns ContinuousStateSpace with NX = DenDeg states.
template<typename Scalar, std::size_t NumDeg, std::size_t DenDeg>
    requires (NumDeg <= DenDeg)
constexpr auto tf2ss(const TransferFunction<Scalar, NumDeg, DenDeg>& tf)
    -> ContinuousStateSpace<Scalar, DenDeg, 1, 1>
{
    static_assert(DenDeg >= 1, "Denominator degree must be at least 1");

    using mat_a = Matrix<Scalar, DenDeg, DenDeg>;
    using mat_b = Matrix<Scalar, DenDeg, 1>;
    using mat_c = Matrix<Scalar, 1, DenDeg>;
    using mat_d = Matrix<Scalar, 1, 1>;

    constexpr auto n = DenDeg;
    constexpr int ni = static_cast<int>(n);

    // Normalize denominator to monic (leading coefficient = 1)
    Scalar a0 = tf.denominator[0];

    std::array<Scalar, n + 1> den{};
    for (std::size_t i = 0; i <= n; ++i)
        den[i] = tf.denominator[i] / a0;

    // Pad numerator to length n+1 (right-aligned, higher degrees get zeros)
    std::array<Scalar, n + 1> num{};
    std::size_t offset = n - NumDeg;
    for (std::size_t i = 0; i <= NumDeg; ++i)
        num[offset + i] = tf.numerator[i] / a0;

    // Controllable canonical form (MATLAB convention)
    mat_a A = mat_a::Zero();

    // Superdiagonal of 1s
    for (int i = 0; i + 1 < ni; ++i)
        A(i, i + 1) = Scalar{1};

    // Last row: [-a_n, -a_{n-1}, ..., -a_1]
    for (int j = 0; j < ni; ++j)
        A(ni - 1, j) = -den[n - static_cast<std::size_t>(j)];

    // B = [0, 0, ..., 1]^T
    mat_b B = mat_b::Zero();
    B(ni - 1, 0) = Scalar{1};

    // C and D depend on whether NumDeg == DenDeg
    mat_c C = mat_c::Zero();
    mat_d D = mat_d::Zero();

    if constexpr (NumDeg == DenDeg) {
        // D = b_0 (leading numerator coefficient, normalized)
        D(0, 0) = num[0];
        // C[j] = num[n-j] - num[0]*den[n-j] for j = 0..n-1
        for (int j = 0; j < ni; ++j)
            C(0, j) = num[n - static_cast<std::size_t>(j)] - num[0] * den[n - static_cast<std::size_t>(j)];
    } else {
        // Strictly proper: D = 0
        D(0, 0) = Scalar{0};
        // C[j] = num[n-j] for j = 0..n-1
        for (int j = 0; j < ni; ++j)
            C(0, j) = num[n - static_cast<std::size_t>(j)];
    }

    return {A, B, C, D};
}

// ss2tf: State-space to transfer function via Leverrier-Faddeev algorithm.
// Computes H(s) = C*(sI - A)^{-1}*B + D for SISO systems.
// Returns TransferFunction<Scalar, NX, NX> (numerator degree = denominator degree = NX).
// Coefficients are highest-degree-first.
template<typename Scalar, std::size_t NX>
auto ss2tf(const ContinuousStateSpace<Scalar, NX, 1, 1>& sys)
    -> TransferFunction<Scalar, NX, NX>
{
    constexpr auto n = NX;
    using mat_nn = Matrix<Scalar, n, n>;

    std::array<Scalar, n + 1> den{};
    std::array<Scalar, n + 1> numer{};

    den[0] = Scalar{1};

    Scalar d_val = sys.D(0, 0);
    numer[0] = d_val;

    mat_nn M = mat_nn::Identity();

    for (std::size_t k = 0; k < n; ++k) {
        // C * M_k * B
        auto CMkB = (sys.C * M * sys.B).eval();
        Scalar cmkb = CMkB(0, 0);

        // A * M_k
        mat_nn AM = (sys.A * M).eval();

        // p_{k+1} = -trace(A * M_k) / (k+1)
        Scalar pk1 = -AM.trace() / static_cast<Scalar>(k + 1);
        den[k + 1] = pk1;

        // numer[k+1] = C*M_k*B + D*p_{k+1}
        numer[k + 1] = cmkb + d_val * pk1;

        // M_{k+1} = A*M_k + p_{k+1}*I
        if (k + 1 < n) {
            M = (AM + pk1 * mat_nn::Identity()).eval();
        }
    }

    return {numer, den};
}

}

#endif
