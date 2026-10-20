#include "ctrlpp/model/conversion.h"

#include "ctrlpp/detail/diagonal_balance.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>

#include <array>
#include <cmath>
#include <limits>
#include <complex>
#include <cstddef>

namespace
{

using Catch::Matchers::WithinAbs;

constexpr double unit_roundoff = std::numeric_limits<double>::epsilon() / 2;

// Higham 2002, Lemma 3.1: k roundings compound to a relative error of at most gamma_k = k u / (1 - k u).
auto gamma(int roundings) -> double
{
    const double ku = roundings * unit_roundoff;
    return ku / (1 - ku);
}

template <int N>
using coefficients = std::array<double, std::size_t(N) + 1>;

template <int N>
auto abs_root_polynomial(const Eigen::Matrix<double, N, 1>& roots) -> coefficients<N>
{
    coefficients<N> c{};
    c[0] = 1;
    for(int i = 0; i < N; ++i)
        for(std::size_t k = std::size_t(N); k >= 1; --k)
            c[k] += roots(i) * c[k - 1];
    return c;
}

template <int N>
struct polynomial_bound
{
    coefficients<N> magnitude;
    coefficients<N> error;
};

// Per coefficient term, the library's factor product rounds twice per real linear factor (product,
// difference) and five times per quadratic factor (a^2 + b^2 and its sum, the product, two additions).
// The bound's own products of positive terms round twice per factor.
constexpr int linear_factor_roundings = 2;
constexpr int quadratic_factor_roundings = 5;
constexpr int bound_roundings_per_factor = 2;

// Golub & Van Loan Sec. 7.5.6: the QR algorithm's eigenvalues are exact for M_bal + E with
// ||E|| <= n u ||M_bal||_F, plus the entrywise formation error of M carried through the same similarity.
// First order, |d lambda_i| <= kappa_i ||E|| with kappa_i = ||y_i|| ||x_i|| / |y_i^H x_i|. With
// rho_i = |lambda_hat_i| + r_i bounding |lambda_i|, every coefficient moves by at most
// e_k(rho + r) - e_k(rho), and the product's rounding adds gamma_m e_k(rho + r) (Higham 2002 Sec. 3.1).
template <int N>
auto characteristic_polynomial_bound(const Eigen::Matrix<double, N, N>& M, const Eigen::Matrix<double, N, N>& formation)
    -> polynomial_bound<N>
{
    Eigen::Matrix<double, N, N> balanced = M;
    Eigen::Matrix<double, N, 1> scaling;
    ctrlpp::detail::balance_diagonal<double, N>(balanced, scaling);
    const double backward = N * unit_roundoff * balanced.norm()
                          + (scaling.cwiseInverse().asDiagonal() * formation * scaling.asDiagonal()).norm();
    const Eigen::EigenSolver<Eigen::Matrix<double, N, N>> right(balanced);
    const Eigen::EigenSolver<Eigen::Matrix<double, N, N>> left(balanced.transpose());
    REQUIRE(right.info() == Eigen::Success);
    REQUIRE(left.info() == Eigen::Success);
    // eigenvectors() returns by value; a column taken straight from it dangles once the statement ends.
    using eigenvectors = typename Eigen::EigenSolver<Eigen::Matrix<double, N, N>>::EigenvectorsType;
    const eigenvectors right_vectors = right.eigenvectors();
    const eigenvectors left_vectors = left.eigenvectors();
    Eigen::Matrix<double, N, 1> reach;
    Eigen::Matrix<double, N, 1> radius;
    int roundings = 0;
    for(int i = 0; i < N; ++i)
    {
        const std::complex<double> lambda = right.eigenvalues()(i);
        Eigen::Index j = 0;
        (left.eigenvalues().array() - lambda).abs().minCoeff(&j);
        const auto x = right_vectors.col(i);
        const auto y = left_vectors.col(j);
        radius(i) = x.norm() * y.norm() / std::abs((y.transpose() * x).value()) * backward;
        reach(i) = std::abs(lambda) + radius(i);
        if(lambda.imag() == 0)
            roundings += linear_factor_roundings;
        else if(lambda.imag() > 0)
            roundings += quadratic_factor_roundings;
    }
    const auto at_reach = abs_root_polynomial<N>(reach);
    const auto beyond = abs_root_polynomial<N>(Eigen::Matrix<double, N, 1>{reach + radius});
    const double evaluation = 2 * gamma(bound_roundings_per_factor * N);
    polynomial_bound<N> out{at_reach, {}};
    for(std::size_t k = 0; k <= std::size_t(N); ++k)
        out.error[k] = (beyond[k] - at_reach[k]) + (gamma(roundings) + evaluation) * beyond[k];
    return out;
}

template <int N>
struct transfer_bound
{
    coefficients<N> numerator;
    coefficients<N> denominator;
};

// num_k = (p1_k - den_k) / alpha + D den_k: the two polynomials' errors propagate through it, and the
// difference, the product and the sum round once each (alpha is a power of two, so the division is exact).
template <int N>
auto transfer_function_bound(const ctrlpp::continuous_state_space<double, std::size_t(N), 1, 1>& sys) -> transfer_bound<N>
{
    using matrix = Eigen::Matrix<double, N, N>;
    const auto den = characteristic_polynomial_bound<N>(sys.A, matrix::Zero());
    const double alpha = ctrlpp::detail::rank_one_scale(sys.A.reshaped().stableNorm(), sys.B.stableNorm(), sys.C.stableNorm());
    const matrix shifted = sys.A - alpha * sys.B * sys.C;
    const matrix formation = gamma(2) * (sys.A.cwiseAbs() + alpha * sys.B.cwiseAbs() * sys.C.cwiseAbs());
    const auto p1 = characteristic_polynomial_bound<N>(shifted, formation);
    const double d = std::abs(sys.D(0, 0));
    transfer_bound<N> out{{}, den.error};
    for(std::size_t k = 1; k <= std::size_t(N); ++k)
    {
        const double p1_reach = p1.magnitude[k] + p1.error[k];
        const double den_reach = den.magnitude[k] + den.error[k];
        out.numerator[k] = (p1.error[k] + den.error[k]) / alpha + d * den.error[k]
                         + gamma(2) * ((p1_reach + den_reach) / alpha + d * den_reach);
    }
    return out;
}

template <std::size_t N>
auto monic_from_real_poles(const std::array<double, N>& poles) -> std::array<double, N + 1>
{
    std::array<double, N + 1> c{};
    c[0] = 1;
    for(std::size_t j = 0; j < N; ++j)
        for(std::size_t k = j + 1; k >= 1; --k)
            c[k] -= poles[j] * c[k - 1];
    return c;
}

template <std::size_t N>
auto ramp_over_poles(const std::array<double, N>& poles) -> ctrlpp::transfer_function<double, N - 1, N>
{
    ctrlpp::transfer_function<double, N - 1, N> tf{};
    for(std::size_t i = 0; i < N; ++i)
        tf.numerator[i] = double(i + 1);
    tf.denominator = monic_from_real_poles(poles);
    return tf;
}

// tf2ss of a monic denominator is exact (negation and copies), so the recovered polynomials are
// compared against the transfer function's own coefficients.
template <std::size_t N>
auto check_round_trip(const ctrlpp::transfer_function<double, N - 1, N>& tf) -> void
{
    auto const realized = ctrlpp::tf2ss(tf);
    REQUIRE(realized.has_value());
    auto const recovered = ctrlpp::ss2tf(*realized);
    REQUIRE(recovered.has_value());
    auto const bound = transfer_function_bound<int(N)>(*realized);

    CHECK(recovered->numerator[0] == realized->D(0, 0));
    for(std::size_t k = 0; k <= N; ++k)
        CHECK_THAT(recovered->denominator[k], WithinAbs(tf.denominator[k], bound.denominator[k]));
    for(std::size_t k = 1; k <= N; ++k)
        CHECK_THAT(recovered->numerator[k], WithinAbs(tf.numerator[k - 1], bound.numerator[k]));
}

}

TEST_CASE("tf2ss converts first-order transfer function")
{
    ctrlpp::transfer_function<double, 0, 1> tf{.numerator = {1.0}, .denominator = {1.0, 2.0}};

    auto const realized = ctrlpp::tf2ss(tf);
    REQUIRE(realized.has_value());
    auto const& ss = *realized;

    CHECK_THAT(ss.A(0, 0), WithinAbs(-2.0, gamma(1) * 2.0));
    CHECK(ss.B(0, 0) == 1.0);
    CHECK_THAT(ss.C(0, 0), WithinAbs(1.0, gamma(1) * 1.0));
    CHECK(ss.D(0, 0) == 0.0);
}

TEST_CASE("tf2ss converts second-order transfer function")
{
    ctrlpp::transfer_function<double, 1, 2> tf{.numerator = {1.0, 1.0}, .denominator = {1.0, 3.0, 2.0}};

    auto const realized = ctrlpp::tf2ss(tf);
    REQUIRE(realized.has_value());
    auto const& ss = *realized;

    CHECK(ss.A(0, 0) == 0.0);
    CHECK(ss.A(0, 1) == 1.0);
    CHECK_THAT(ss.A(1, 0), WithinAbs(-2.0, gamma(1) * 2.0));
    CHECK_THAT(ss.A(1, 1), WithinAbs(-3.0, gamma(1) * 3.0));

    CHECK(ss.B(0, 0) == 0.0);
    CHECK(ss.B(1, 0) == 1.0);

    CHECK_THAT(ss.C(0, 0), WithinAbs(1.0, gamma(1) * 1.0));
    CHECK_THAT(ss.C(0, 1), WithinAbs(1.0, gamma(1) * 1.0));

    CHECK(ss.D(0, 0) == 0.0);
}

TEST_CASE("ss2tf recovers transfer function coefficients")
{
    using SS = ctrlpp::continuous_state_space<double, 1, 1, 1>;
    SS sys{};
    sys.A(0, 0) = -2.0;
    sys.B(0, 0) = 1.0;
    sys.C(0, 0) = 1.0;
    sys.D(0, 0) = 0.0;

    auto const tf = ctrlpp::ss2tf(sys);
    REQUIRE(tf.has_value());
    auto const bound = transfer_function_bound<1>(sys);

    CHECK(tf->numerator[0] == 0.0);
    CHECK_THAT(tf->numerator[1], WithinAbs(1.0, bound.numerator[1]));
    CHECK(tf->denominator[0] == 1.0);
    CHECK_THAT(tf->denominator[1], WithinAbs(2.0, bound.denominator[1]));
}

TEST_CASE("tf2ss -> ss2tf roundtrip preserves coefficients")
{
    check_round_trip<2>({.numerator = {1.0, 1.0}, .denominator = {1.0, 3.0, 2.0}});
}

TEST_CASE("tf2ss handles proper transfer function with equal degrees")
{
    ctrlpp::transfer_function<double, 1, 1> tf{.numerator = {2.0, 5.0}, .denominator = {1.0, 3.0}};

    auto const realized = ctrlpp::tf2ss(tf);
    REQUIRE(realized.has_value());
    auto const& ss = *realized;

    CHECK_THAT(ss.A(0, 0), WithinAbs(-3.0, gamma(1) * 3.0));
    CHECK(ss.B(0, 0) == 1.0);
    CHECK_THAT(ss.C(0, 0), WithinAbs(-1.0, gamma(4) * (5.0 + 2.0 * 3.0)));
    CHECK_THAT(ss.D(0, 0), WithinAbs(2.0, gamma(1) * 2.0));
}

TEST_CASE("tf2ss refuses a zero leading denominator coefficient")
{
    ctrlpp::transfer_function<double, 0, 2> tf{.numerator = {1.0}, .denominator = {0.0, 1.0, 2.0}};

    auto const realized = ctrlpp::tf2ss(tf);

    REQUIRE_FALSE(realized.has_value());
    CHECK(realized.error() == ctrlpp::conversion_error::zero_leading_coefficient);
}

TEST_CASE("tf2ss refuses non-finite coefficients")
{
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    constexpr double inf = std::numeric_limits<double>::infinity();
    ctrlpp::transfer_function<double, 1, 2> nan_numerator{.numerator = {1.0, nan}, .denominator = {1.0, 3.0, 2.0}};
    ctrlpp::transfer_function<double, 1, 2> inf_denominator{.numerator = {1.0, 1.0}, .denominator = {1.0, inf, 2.0}};

    auto const from_nan = ctrlpp::tf2ss(nan_numerator);
    auto const from_inf = ctrlpp::tf2ss(inf_denominator);

    REQUIRE_FALSE(from_nan.has_value());
    CHECK(from_nan.error() == ctrlpp::conversion_error::non_finite_input);
    REQUIRE_FALSE(from_inf.has_value());
    CHECK(from_inf.error() == ctrlpp::conversion_error::non_finite_input);
}

TEST_CASE("tf2ss refuses a leading coefficient whose quotients overflow")
{
    ctrlpp::transfer_function<double, 0, 2> tf{.numerator = {1.0}, .denominator = {1e-300, 1e10, 1.0}};

    auto const realized = ctrlpp::tf2ss(tf);

    REQUIRE_FALSE(realized.has_value());
    CHECK(realized.error() == ctrlpp::conversion_error::unrepresentable_result);
}

TEST_CASE("ss2tf recovers eight poles spread over four decades")
{
    check_round_trip<8>(ramp_over_poles<8>({-0.1, -0.2, -0.5, -1.0, -2.0, -5.0, -10.0, -1000.0}));
}

TEST_CASE("ss2tf recovers six poles spread over three decades")
{
    check_round_trip<6>(ramp_over_poles<6>({-0.1, -0.2, -1.0, -5.0, -20.0, -100.0}));
}

TEST_CASE("ss2tf recovers eight clustered poles")
{
    check_round_trip<8>(ramp_over_poles<8>({-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0}));
}

TEST_CASE("ss2tf recovers four poles spread over three decades")
{
    check_round_trip<4>(ramp_over_poles<4>({-0.1, -1.0, -10.0, -100.0}));
}

TEST_CASE("ss2tf recovers two lightly damped conjugate pairs")
{
    // (s^2 + 0.2 s + 4)(s^2 + 0.5 s + 100)
    check_round_trip<4>({.numerator = {1.0, 2.0, 3.0, 4.0}, .denominator = {1.0, 0.7, 104.1, 22.0, 400.0}});
}

TEST_CASE("ss2tf refuses a non-finite state-space entry")
{
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    auto const realized = ctrlpp::tf2ss(ctrlpp::transfer_function<double, 1, 2>{.numerator = {1.0, 1.0}, .denominator = {1.0, 3.0, 2.0}});
    REQUIRE(realized.has_value());

    for(int entry = 0; entry < 4; ++entry)
    {
        auto sys = *realized;
        if(entry == 0)
            sys.A(1, 0) = nan;
        else if(entry == 1)
            sys.B(0, 0) = nan;
        else if(entry == 2)
            sys.C(0, 1) = nan;
        else
            sys.D(0, 0) = nan;

        auto const tf = ctrlpp::ss2tf(sys);

        REQUIRE_FALSE(tf.has_value());
        CHECK(tf.error() == ctrlpp::conversion_error::non_finite_input);
    }
}

TEST_CASE("ss2tf with a zero input or output map returns D times the denominator")
{
    auto const realized = ctrlpp::tf2ss(ctrlpp::transfer_function<double, 2, 2>{.numerator = {3.0, 1.0, 1.0}, .denominator = {1.0, 3.0, 2.0}});
    REQUIRE(realized.has_value());

    for(bool zero_input : {true, false})
    {
        auto sys = *realized;
        if(zero_input)
            sys.B.setZero();
        else
            sys.C.setZero();

        auto const tf = ctrlpp::ss2tf(sys);

        REQUIRE(tf.has_value());
        for(std::size_t k = 0; k <= 2; ++k)
            CHECK(tf->numerator[k] == sys.D(0, 0) * tf->denominator[k]);
    }
}
