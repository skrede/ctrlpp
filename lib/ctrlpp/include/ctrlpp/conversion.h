#ifndef HPP_GUARD_CTRLPP_CONVERSION_H
#define HPP_GUARD_CTRLPP_CONVERSION_H

#include "ctrlpp/state_space.h"
#include "ctrlpp/transfer_function.h"

#include <array>
#include <cstddef>

namespace ctrlpp {

namespace detail {

template<typename Scalar, std::size_t N, typename Policy>
concept ElementAccessPolicy = LinalgPolicy<Policy> && requires(
    typename Policy::template matrix_type<Scalar, N, N> m,
    Scalar val
) {
    { Policy::get_element(m, std::size_t{}, std::size_t{}) } -> std::convertible_to<Scalar>;
    Policy::set_element(m, std::size_t{}, std::size_t{}, val);
};

template<typename Scalar, std::size_t N, typename Policy>
concept TracePolicy = requires(
    typename Policy::template matrix_type<Scalar, N, N> m
) {
    { Policy::trace(m) } -> std::convertible_to<Scalar>;
};

}

// tf2ss: Transfer function to controllable canonical form state-space.
// H(s) = num(s) / den(s), coefficients highest-degree-first (MATLAB convention).
// Requires NumDeg <= DenDeg (proper transfer function).
// Returns ContinuousStateSpace with NX = DenDeg states.
template<typename Scalar, std::size_t NumDeg, std::size_t DenDeg, LinalgPolicy Policy>
    requires (NumDeg <= DenDeg) &&
             detail::ElementAccessPolicy<Scalar, DenDeg, Policy>
constexpr auto tf2ss(const TransferFunction<Scalar, NumDeg, DenDeg, Policy>& tf)
    -> ContinuousStateSpace<Scalar, DenDeg, 1, 1, Policy>
{
    static_assert(DenDeg >= 1, "Denominator degree must be at least 1");

    using mat_a = typename Policy::template matrix_type<Scalar, DenDeg, DenDeg>;
    using mat_b = typename Policy::template matrix_type<Scalar, DenDeg, 1>;
    using mat_c = typename Policy::template matrix_type<Scalar, 1, DenDeg>;
    using mat_d = typename Policy::template matrix_type<Scalar, 1, 1>;

    constexpr auto n = DenDeg;

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

    // Controllable canonical form (MATLAB convention):
    // A = [[0, 1, 0, ..., 0    ],
    //      [0, 0, 1, ..., 0    ],
    //      [...                 ],
    //      [0, 0, 0, ..., 1    ],
    //      [-a_n, -a_{n-1}, ..., -a_1]]
    // where den = [1, a_1, a_2, ..., a_n] (monic, highest-degree-first)
    mat_a A{};
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            Policy::set_element(A, i, j, Scalar{0});

    // Superdiagonal of 1s
    for (std::size_t i = 0; i + 1 < n; ++i)
        Policy::set_element(A, i, i + 1, Scalar{1});

    // Last row: [-a_n, -a_{n-1}, ..., -a_1]
    for (std::size_t j = 0; j < n; ++j)
        Policy::set_element(A, n - 1, j, -den[n - j]);

    // B = [0, 0, ..., 1]^T
    mat_b B{};
    for (std::size_t i = 0; i < n; ++i)
        Policy::set_element(B, i, 0, Scalar{0});
    Policy::set_element(B, n - 1, 0, Scalar{1});

    // C and D depend on whether NumDeg == DenDeg
    mat_c C{};
    mat_d D{};

    if constexpr (NumDeg == DenDeg) {
        // D = b_0 (leading numerator coefficient, normalized)
        Policy::set_element(D, 0, 0, num[0]);
        // C[j] = num[n-j] - num[0]*den[n-j] for j = 0..n-1
        for (std::size_t j = 0; j < n; ++j)
            Policy::set_element(C, 0, j, num[n - j] - num[0] * den[n - j]);
    } else {
        // Strictly proper: D = 0
        Policy::set_element(D, 0, 0, Scalar{0});
        // C[j] = num[n-j] for j = 0..n-1
        for (std::size_t j = 0; j < n; ++j)
            Policy::set_element(C, 0, j, num[n - j]);
    }

    return {A, B, C, D};
}

// ss2tf: State-space to transfer function via Leverrier-Faddeev algorithm.
// Computes H(s) = C*(sI - A)^{-1}*B + D for SISO systems.
// Returns TransferFunction<Scalar, NX, NX, Policy> (numerator degree = denominator degree = NX).
// Coefficients are highest-degree-first.
template<typename Scalar, std::size_t NX, LinalgPolicy Policy>
    requires detail::ElementAccessPolicy<Scalar, NX, Policy> &&
             detail::TracePolicy<Scalar, NX, Policy>
auto ss2tf(const ContinuousStateSpace<Scalar, NX, 1, 1, Policy>& sys)
    -> TransferFunction<Scalar, NX, NX, Policy>
{
    constexpr auto n = NX;
    using mat_nn = typename Policy::template matrix_type<Scalar, n, n>;

    // Leverrier-Faddeev algorithm:
    // adj(sI - A) = M_0*s^{n-1} + M_1*s^{n-2} + ... + M_{n-1}
    // det(sI - A) = s^n + p_1*s^{n-1} + ... + p_n
    //
    // Recursion: M_0 = I, then for k = 0..n-1:
    //   p_{k+1} = -trace(A * M_k) / (k+1)
    //   M_{k+1} = A*M_k + p_{k+1}*I  (only needed for k+1 < n)
    //
    // H(s) = (C*adj(sI-A)*B + D*det(sI-A)) / det(sI-A)
    //
    // Numerator coefficients (highest-degree-first):
    //   numer[0]   = D                           (coeff of s^n)
    //   numer[k+1] = C*M_k*B + D*p_{k+1}        (coeff of s^{n-k-1}, for k=0..n-1)

    std::array<Scalar, n + 1> den{};
    std::array<Scalar, n + 1> numer{};

    den[0] = Scalar{1};

    Scalar d_val = Policy::get_element(sys.D, std::size_t{0}, std::size_t{0});
    numer[0] = d_val;

    mat_nn M = Policy::template identity<Scalar, n>();

    for (std::size_t k = 0; k < n; ++k) {
        // C * M_k * B
        auto CMkB = Policy::multiply(sys.C, Policy::multiply(M, sys.B));
        Scalar cmkb = Policy::get_element(CMkB, std::size_t{0}, std::size_t{0});

        // A * M_k
        mat_nn AM = Policy::multiply(sys.A, M);

        // p_{k+1} = -trace(A * M_k) / (k+1)
        Scalar pk1 = -Policy::trace(AM) / static_cast<Scalar>(k + 1);
        den[k + 1] = pk1;

        // numer[k+1] = C*M_k*B + D*p_{k+1}
        numer[k + 1] = cmkb + d_val * pk1;

        // M_{k+1} = A*M_k + p_{k+1}*I
        if (k + 1 < n) {
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j) {
                    Scalar val = Policy::get_element(AM, i, j);
                    if (i == j)
                        val += pk1;
                    Policy::set_element(M, i, j, val);
                }
        }
    }

    return {numer, den};
}

}

#endif
