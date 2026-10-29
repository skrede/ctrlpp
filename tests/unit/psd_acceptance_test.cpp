#include "ctrlpp/detail/psd_acceptance.h"
#include "ctrlpp/detail/power_of_two_scale.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/QR>
#include <Eigen/Core>

#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <cstdint>
#include <algorithm>

namespace
{

template<typename Scalar, int N>
using square = Eigen::Matrix<Scalar, N, N>;

constexpr int draws = 1000;

constexpr std::array<int, 3> member_exponents{0, -1000, 1000};

// Raw engine output mapped by explicit arithmetic, because the standard distributions differ between standard libraries.
auto unit_interval(std::mt19937_64 &engine) -> double
{
    return static_cast<double>(engine() >> 11) * 0x1p-53;
}

auto small_integer(std::mt19937_64 &engine) -> double
{
    return static_cast<double>(engine() % 9) - 4.0;
}

auto random_sign(std::mt19937_64 &engine) -> double
{
    return (engine() & 1U) != 0 ? 1.0 : -1.0;
}

// Each entry is fl(v_i v_j), one rounding of the exact rank-one v v'. |v_i| >= 1/8 and the grading below keep every
// member and its 2^-1000 multiple in the normal range, so a power-of-two multiple is still one rounding per entry.
template<typename Scalar, int N>
auto rounded_rank_one(std::mt19937_64 &engine) -> square<Scalar, N>
{
    Eigen::Matrix<Scalar, N, 1> v;
    for(int i = 0; i < N; ++i)
        v(i) = static_cast<Scalar>(random_sign(engine) * (0.125 + 0.875 * unit_interval(engine)));
    square<Scalar, N> P;
    for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j)
            P(i, j) = v(i) * v(j);
    return P;
}

// A sum of outer products of small integer vectors is exact, so fl(c G) with c = fl(1/10) rounds each entry of the
// semi-definite c G once.
template<int N>
auto rounded_rank(std::mt19937_64 &engine, int rank) -> square<double, N>
{
    square<double, N> G = square<double, N>::Zero();
    for(int k = 0; k < rank; ++k)
    {
        Eigen::Matrix<double, N, 1> w;
        for(int i = 0; i < N; ++i)
            w(i) = small_integer(engine);
        G += w * w.transpose();
    }
    return 0.1 * G;
}

template<int N>
auto rounded_member(std::mt19937_64 &engine, int rank) -> square<double, N>
{
    return rank == 1 ? rounded_rank_one<double, N>(engine) : rounded_rank<N>(engine, rank);
}

template<int N>
auto graded(const square<double, N> &P, std::mt19937_64 &engine) -> square<double, N>
{
    Eigen::Matrix<double, N, 1> g;
    for(int i = 0; i < N; ++i)
        g(i) = std::ldexp(1.0, -static_cast<int>(engine() % 7));
    return (P.array() * (g * g.transpose()).array()).matrix();
}

template<int N>
auto aligned(square<double, N> P, std::mt19937_64 &engine) -> square<double, N>
{
    const std::uint64_t zeroed = 1 + engine() % ((std::uint64_t{1} << N) - 2);
    for(int i = 0; i < N; ++i)
        if(((zeroed >> i) & 1U) != 0)
        {
            P.row(i).setZero();
            P.col(i).setZero();
        }
    return P;
}

template<int N>
auto shaped(const square<double, N> &P, int shape, std::mt19937_64 &engine) -> square<double, N>
{
    if(shape == 1)
        return graded<N>(P, engine);
    if(shape == 2)
        return aligned<N>(P, engine);
    return P;
}

// H / 2 with H the 4 x 4 Hadamard matrix is orthogonal and exact; other sizes take the Q of a Householder QR.
template<int N>
auto rotation(std::mt19937_64 &engine) -> square<double, N>
{
    if constexpr(N == 4)
    {
        square<double, 4> H;
        H << 1, 1, 1, 1, 1, -1, 1, -1, 1, 1, -1, -1, 1, -1, -1, 1;
        return H / 2.0;
    }
    else
    {
        square<double, N> M;
        for(int i = 0; i < N; ++i)
            for(int j = 0; j < N; ++j)
                M(i, j) = 2.0 * unit_interval(engine) - 1.0;
        return Eigen::HouseholderQR<square<double, N>>(M).householderQ();
    }
}

template<int N>
auto rotated(const square<double, N> &P, std::mt19937_64 &engine) -> square<double, N>
{
    const square<double, N> Q = rotation<N>(engine);
    const square<double, N> R = Q * P * Q.transpose();
    return (R + R.transpose()) / 2.0;
}

template<int N>
auto far_indefinite(double delta, int position) -> square<double, N>
{
    square<double, N> P   = square<double, N>::Identity();
    P(position, position) = -delta;
    return P;
}

// With A = diag(a_i) and B = Q = R = I each mode's Riccati roots are a_i +- sqrt(a_i^2 + 1); the minus root on one
// mode gives a mixed, indefinite solution.
template<int N>
auto mixed_riccati_solution(std::mt19937_64 &engine, int swapped) -> square<double, N>
{
    square<double, N> P = square<double, N>::Zero();
    for(int i = 0; i < N; ++i)
    {
        const double a    = 4.0 * unit_interval(engine) - 2.0;
        const double root = std::sqrt(a * a + 1.0);
        P(i, i)           = i == swapped ? a - root : a + root;
    }
    return P;
}

// The remaining diagonal is `rest` times the block's own magnitude, so a small coupling stays far beyond rounding.
template<int N>
auto coupled(double diagonal, double coupling, int row, int col, double rest) -> square<double, N>
{
    square<double, N> P = rest * std::max(diagonal, coupling) * square<double, N>::Identity();
    P(row, row)         = diagonal;
    P(col, col)         = 0.0;
    P(row, col)         = coupling;
    P(col, row)         = coupling;
    return P;
}

template<int N>
auto check_rounded_family(std::mt19937_64 &engine, int rank) -> void
{
    constexpr std::array<const char *, 3> shapes{"plain", "graded", "coordinate-aligned"};
    for(int shape = 0; shape < 3; ++shape)
    {
        std::array<int, 3> refused{};
        for(int draw = 0; draw < draws; ++draw)
        {
            const square<double, N> member = shaped<N>(rounded_member<N>(engine, rank), shape, engine);
            for(std::size_t e = 0; e < member_exponents.size(); ++e)
                if(!ctrlpp::detail::psd_within_rounding<double, N>(std::ldexp(1.0, member_exponents[e]) * member))
                    ++refused[e];
        }
        for(std::size_t e = 0; e < member_exponents.size(); ++e)
        {
            CAPTURE(N, rank, shapes[static_cast<std::size_t>(shape)], member_exponents[e]);
            CHECK(refused[e] == 0);
        }
    }
}

template<int N>
auto check_rounded_ranks(std::uint64_t seed) -> void
{
    std::mt19937_64 engine(seed);
    for(int rank = 1; rank <= N; ++rank)
        check_rounded_family<N>(engine, rank);
}

template<int N>
auto check_zero_block_coupling() -> void
{
    constexpr std::array<double, 5> zero_diagonal{std::numeric_limits<double>::denorm_min(), 0x1p-1000, 1e-300, 1.0, 1e300};
    constexpr std::array<double, 3> unit_diagonal{1e-6, 1e-3, 0.5};
    int drawn   = 0;
    int refused = 0;
    for(int row = 0; row < N; ++row)
        for(int col = 0; col < N; ++col)
            for(const double rest : {0.0, 1.0})
            {
                if(row == col || (N == 2 && rest != 0.0))
                    continue;
                for(const double c : zero_diagonal)
                    refused += ctrlpp::detail::psd_within_rounding<double, N>(coupled<N>(0.0, c, row, col, rest)) ? 0 : 1;
                for(const double c : unit_diagonal)
                    refused += ctrlpp::detail::psd_within_rounding<double, N>(coupled<N>(1.0, c, row, col, rest)) ? 0 : 1;
                drawn += static_cast<int>(zero_diagonal.size() + unit_diagonal.size());
            }
    CAPTURE(N);
    CHECK(refused == drawn);
}

template<int N>
auto check_mixed_solutions(std::uint64_t seed) -> void
{
    std::mt19937_64 engine(seed);
    std::array<int, 2> refused{};
    for(int draw = 0; draw < draws; ++draw)
    {
        const square<double, N> P = mixed_riccati_solution<N>(engine, draw % N);
        refused[0] += ctrlpp::detail::psd_within_rounding<double, N>(P) ? 0 : 1;
        refused[1] += ctrlpp::detail::psd_within_rounding<double, N>(rotated<N>(P, engine)) ? 0 : 1;
    }
    CAPTURE(N);
    CHECK(refused[0] == draws);
    CHECK(refused[1] == draws);
}

template<int N>
auto check_far_indefinite(std::uint64_t seed) -> void
{
    std::mt19937_64 engine(seed);
    for(const double delta : {1e-8, 1e-4})
    {
        std::array<int, 2> refused{};
        for(int draw = 0; draw < draws; ++draw)
        {
            const square<double, N> P = far_indefinite<N>(delta, draw % N);
            refused[0] += ctrlpp::detail::psd_within_rounding<double, N>(P) ? 0 : 1;
            refused[1] += ctrlpp::detail::psd_within_rounding<double, N>(rotated<N>(P, engine)) ? 0 : 1;
        }
        CAPTURE(N, delta);
        CHECK(refused[0] == draws);
        CHECK(refused[1] == draws);
    }
}

template<int N>
auto check_non_finite_and_zero() -> void
{
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    constexpr double inf = std::numeric_limits<double>::infinity();
    int drawn            = 0;
    int refused          = 0;
    for(const double bad : {nan, inf, -inf})
    {
        square<double, N> diagonal = square<double, N>::Identity();
        diagonal(N - 1, N - 1)     = bad;
        square<double, N> coupling = square<double, N>::Identity();
        coupling(0, N - 1)         = bad;
        coupling(N - 1, 0)         = bad;
        refused += ctrlpp::detail::psd_within_rounding<double, N>(diagonal) ? 0 : 1;
        refused += ctrlpp::detail::psd_within_rounding<double, N>(coupling) ? 0 : 1;
        drawn += 2;
    }
    CAPTURE(N);
    CHECK(refused == drawn);
    CHECK(ctrlpp::detail::psd_within_rounding<double, N>(square<double, N>::Zero()));
}

template<int N>
auto check_scaling_invariance(std::uint64_t seed) -> void
{
    std::mt19937_64 engine(seed);
    int changed = 0;
    for(int draw = 0; draw < draws; ++draw)
    {
        const int rank            = 1 + draw % N;
        const square<double, N> P = draw % 2 == 0 ? rounded_member<N>(engine, rank) : rotated<N>(far_indefinite<N>(1e-8, draw % N), engine);
        const bool verdict        = ctrlpp::detail::psd_within_rounding<double, N>(P);
        for(const int k : {-1000, -1, 1, 1000})
            changed += ctrlpp::detail::psd_within_rounding<double, N>(std::ldexp(1.0, k) * P) == verdict ? 0 : 1;
    }
    CAPTURE(N);
    CHECK(changed == 0);
}

// The real-valued bound in the theorem's own constants, which the evaluated shift must exceed:
// tau > [u f / (1 - u) + (u + c_n (1 + u)) d] / [1 - u - c_n (1 + u)], c_n = n gamma_{n+1} / (1 - gamma_{n+1}).
template<typename Wide, int N>
auto shift_bound(Wide u, Wide frobenius, Wide largest_diagonal) -> Wide
{
    const Wide n          = static_cast<Wide>(N);
    const Wide gamma      = (n + 1) * u / (1 - (n + 1) * u);
    const Wide completion = n * gamma / (1 - gamma);
    const Wide grown      = completion * (1 + u);
    return (u * frobenius / (1 - u) + (u + grown) * largest_diagonal) / (1 - u - grown);
}

// The zero matrix is accepted before any shift is formed, and its bound is zero, so it is not posed here.
template<typename Scalar, typename Wide, int N>
auto shift_exceeds_bound(const square<Scalar, N> &P) -> bool
{
    if(P.isZero(Scalar{0}))
        return true;
    const square<Scalar, N> scaled = P / ctrlpp::detail::power_of_two_ceiling(P.cwiseAbs().maxCoeff());
    const Scalar largest_diagonal  = std::max(Scalar{0}, scaled.diagonal().maxCoeff());
    const Scalar shift             = ctrlpp::detail::cholesky_completion_shift<Scalar, N>(scaled.norm(), largest_diagonal);
    const Wide u                   = static_cast<Wide>(std::numeric_limits<Scalar>::epsilon()) / 2;
    const Wide frobenius           = scaled.template cast<Wide>().norm();
    return static_cast<Wide>(shift) > shift_bound<Wide, N>(u, frobenius, static_cast<Wide>(largest_diagonal));
}

template<int N>
auto check_double_shift_bound(std::uint64_t seed) -> void
{
    std::mt19937_64 engine(seed);
    int below = 0;
    for(int draw = 0; draw < draws; ++draw)
    {
        const square<double, N> rounded    = shaped<N>(rounded_member<N>(engine, 1 + draw % N), draw % 3, engine);
        const square<double, N> indefinite = rotated<N>(far_indefinite<N>(1e-4, draw % N), engine);
        below += shift_exceeds_bound<double, long double, N>(rounded) ? 0 : 1;
        below += shift_exceeds_bound<double, long double, N>(indefinite) ? 0 : 1;
    }
    CAPTURE(N);
    CHECK(below == 0);
}

template<int N>
auto check_float_shift_bound(std::uint64_t seed) -> void
{
    std::mt19937_64 engine(seed);
    int below = 0;
    for(int draw = 0; draw < draws; ++draw)
        below += shift_exceeds_bound<float, double, N>(rounded_rank_one<float, N>(engine)) ? 0 : 1;
    CAPTURE(N);
    CHECK(below == 0);
}

}

TEST_CASE("Every matrix within one rounding per entry of a positive semi-definite matrix is accepted", "[riccati][definiteness]")
{
    check_rounded_ranks<2>(2);
    check_rounded_ranks<3>(3);
    check_rounded_ranks<4>(4);
    check_rounded_ranks<8>(8);
}

TEST_CASE("A coupling into a zero diagonal block is refused at every scale", "[riccati][definiteness]")
{
    check_zero_block_coupling<2>();
    check_zero_block_coupling<3>();
    check_zero_block_coupling<8>();
}

TEST_CASE("A mixed Riccati solution is refused whether rotated or not", "[riccati][definiteness]")
{
    check_mixed_solutions<2>(12);
    check_mixed_solutions<3>(13);
    check_mixed_solutions<4>(14);
    check_mixed_solutions<8>(18);
}

TEST_CASE("A matrix indefinite far beyond rounding is refused whether rotated or not", "[riccati][definiteness]")
{
    check_far_indefinite<2>(22);
    check_far_indefinite<3>(23);
    check_far_indefinite<4>(24);
    check_far_indefinite<8>(28);
}

TEST_CASE("A non-finite matrix is refused and the zero matrix is accepted", "[riccati][definiteness]")
{
    check_non_finite_and_zero<2>();
    check_non_finite_and_zero<3>();
    check_non_finite_and_zero<4>();
    check_non_finite_and_zero<8>();
}

TEST_CASE("The verdict is unchanged by an exact power-of-two scaling", "[riccati][definiteness]")
{
    check_scaling_invariance<2>(32);
    check_scaling_invariance<3>(33);
    check_scaling_invariance<4>(34);
    check_scaling_invariance<8>(38);
}

TEST_CASE("Single precision accepts a rounded rank-one matrix and refuses a zero-block coupling", "[riccati][definiteness]")
{
    std::mt19937_64 engine(44);
    int refused = 0;
    for(int draw = 0; draw < draws; ++draw)
        refused += ctrlpp::detail::psd_within_rounding<float, 4>(rounded_rank_one<float, 4>(engine)) ? 0 : 1;
    CHECK(refused == 0);

    constexpr std::array<float, 4> couplings{std::numeric_limits<float>::denorm_min(), 0x1p-100F, 1.0F, 0x1p100F};
    int accepted = 0;
    for(const float c : couplings)
    {
        square<float, 2> P;
        P << 0.0F, c, c, 0.0F;
        accepted += ctrlpp::detail::psd_within_rounding<float, 2>(P) ? 1 : 0;
    }
    CHECK(accepted == 0);
}

TEST_CASE("The evaluated shift lies strictly above its real-valued bound", "[riccati][definiteness]")
{
    check_float_shift_bound<2>(62);
    check_float_shift_bound<3>(63);
    check_float_shift_bound<4>(64);
    check_float_shift_bound<8>(68);

    // SKIP needs exceptions, which the default test tree compiles out.
    if(std::numeric_limits<long double>::digits <= std::numeric_limits<double>::digits)
    {
        WARN("long double is no wider than double here, so the double shift is not checked against its bound");
        return;
    }
    check_double_shift_bound<2>(52);
    check_double_shift_bound<3>(53);
    check_double_shift_bound<4>(54);
    check_double_shift_bound<8>(58);
}
