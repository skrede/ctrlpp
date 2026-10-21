// What the oracles in this file decide.
//
// At one output the five Kalman-family filters solve their gain by the pivot
// of the 1 x 1 innovation covariance instead of factoring it. The claim is that
// this is bitwise the column-pivoting QR solve the filters use at every other
// output dimension, so every case compares bits against that solve:
//
//  * The solve itself, on every right-hand form and column count the filters
//    pass, over random operands from the subnormal range to overflow and over
//    the cross product of the special values (signed zeros, infinities, NaNs,
//    subnormals, the extremes). A NaN result is matched as NaN, because which
//    input NaN an operation propagates is unspecified.
//  * Each filter's first update at one output, built so that its innovation
//    covariance and right-hand side are exact (or, for manifold_ukf, recomputed
//    here from the same sigma points), with an innovation of exactly one so the
//    corrected state reads the gain without rounding. Each case built from
//    exact data also checks that its data tells division from multiplication by
//    the reciprocal, so a path that picked the wrong one would fail.
//
// The reference QR is dynamic-size: the triangular kernel Eigen selects depends
// only on the right-hand side's compile-time column count, which is kept.

#include "hardening_helpers.h"

#include "ctrlpp/estimation/ekf.h"
#include "ctrlpp/estimation/ukf.h"
#include "ctrlpp/estimation/mekf.h"
#include "ctrlpp/estimation/kalman.h"
#include "ctrlpp/estimation/manifold_ukf.h"

#include "ctrlpp/estimation/sigma_points/so3_sigma_points.h"
#include "ctrlpp/estimation/sigma_points/julier_sigma_points.h"

#include "ctrlpp/detail/covariance_ops.h"

#include "ctrlpp/lie/so3.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <cstdint>
#include <cstring>

namespace
{

using ctrlpp::Matrix;
using ctrlpp::Vector;

template <typename Scalar, int N>
using row = Eigen::Matrix<Scalar, 1, N>;

template <typename A, typename B>
auto same_bits(const A& a, const B& b) -> bool
{
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), sizeof(typename A::Scalar) * static_cast<std::size_t>(a.size())) == 0;
}

// Two NaN operands of opposite sign give a NaN whose sign is whichever operand
// the instruction propagates, and the operand order of a commutative product is
// the compiler's to choose (IEEE 754-2019, 6.2.3 and 6.3), so it differs even
// between Eigen's own packet and scalar paths. Every other result is compared
// bitwise.
template <typename Scalar, int N>
auto same_solution(const row<Scalar, N>& a, const row<Scalar, N>& b) -> bool
{
    for(int j = 0; j < N; ++j)
    {
        const bool both_nan = std::isnan(a(j)) && std::isnan(b(j));
        if(!both_nan && std::memcmp(&a(j), &b(j), sizeof(Scalar)) != 0)
            return false;
    }
    return true;
}

template <typename Scalar, int N>
auto qr_solve(Scalar s, const row<Scalar, N>& b) -> row<Scalar, N>
{
    using rhs_t = Eigen::Matrix<Scalar, Eigen::Dynamic, N, N == 1 ? Eigen::ColMajor : Eigen::RowMajor>;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> S(1, 1);
    S(0, 0) = s;
    const rhs_t B = b;
    const rhs_t K = S.transpose().colPivHouseholderQr().solve(B);
    return K;
}

template <typename Scalar>
auto specials() -> std::array<Scalar, 16>
{
    using limits = std::numeric_limits<Scalar>;
    return {Scalar{0},         -Scalar{0},         limits::infinity(), -limits::infinity(), limits::quiet_NaN(), -limits::quiet_NaN(),
            limits::denorm_min(), -limits::denorm_min(), limits::min(),   -limits::min(),      Scalar{1},           -Scalar{1},
            Scalar{3},         -Scalar{3},         limits::max(),      -limits::max()};
}

template <typename Scalar>
auto random_operand(std::mt19937_64& rng) -> Scalar
{
    using limits = std::numeric_limits<Scalar>;
    std::uniform_int_distribution<int> exponent(limits::min_exponent - limits::digits, limits::max_exponent);
    std::uniform_real_distribution<Scalar> mantissa(Scalar{1}, Scalar{2});
    std::bernoulli_distribution negative(0.5);
    const Scalar magnitude = std::ldexp(mantissa(rng), exponent(rng));
    return negative(rng) ? -magnitude : magnitude;
}

template <typename Scalar, int N>
auto every_form_matches(Scalar s, const row<Scalar, N>& b) -> bool
{
    const Eigen::Matrix<Scalar, N, 1> column = b.transpose();
    const Eigen::Matrix<Scalar, N, N> identity = Eigen::Matrix<Scalar, N, N>::Identity();
    const row<Scalar, N> product = b * identity;
    const row<Scalar, N> plain_gain = ctrlpp::detail::one_output_gain_transpose(s, b);
    const row<Scalar, N> transpose_gain = ctrlpp::detail::one_output_gain_transpose(s, column.transpose());
    const row<Scalar, N> product_gain = ctrlpp::detail::one_output_gain_transpose(s, b * identity);
    const row<Scalar, N> reference = qr_solve(s, b);
    return same_solution(plain_gain, reference) && same_solution(transpose_gain, reference) && same_solution(product_gain, qr_solve(s, product));
}

template <typename Scalar, int N>
auto mismatches(std::uint64_t seed) -> int
{
    std::mt19937_64 rng(seed);
    int count = 0;
    for(int pair = 0; pair < 4096; ++pair)
    {
        const Scalar s = random_operand<Scalar>(rng);
        row<Scalar, N> b;
        for(int j = 0; j < N; ++j)
            b(j) = random_operand<Scalar>(rng);
        count += every_form_matches(s, b) ? 0 : 1;
    }
    for(const Scalar s : specials<Scalar>())
        for(const Scalar value : specials<Scalar>())
        {
            row<Scalar, N> b = row<Scalar, N>::Constant(value);
            b(N - 1) = N == 1 ? value : Scalar{3};
            count += every_form_matches(s, b) ? 0 : 1;
        }
    return count;
}

struct hold_dynamics
{
    auto operator()(const Vector<double, 2>& x, const Vector<double, 1>&) const -> Vector<double, 2> { return x; }

    auto jacobian_x(const Vector<double, 2>&, const Vector<double, 1>&) const -> Matrix<double, 2, 2> { return Matrix<double, 2, 2>::Identity(); }

    auto jacobian_u(const Vector<double, 2>&, const Vector<double, 1>&) const -> Matrix<double, 2, 1> { return Matrix<double, 2, 1>::Zero(); }
};

template <std::size_t NX>
struct first_state
{
    auto operator()(const Vector<double, NX>& x) const -> Vector<double, 1> { return x.template head<1>(); }

    auto jacobian(const Vector<double, NX>&) const -> Matrix<double, 1, NX>
    {
        Matrix<double, 1, NX> H = Matrix<double, 1, NX>::Zero();
        H(0, 0) = 1.0;
        return H;
    }
};

template <std::size_t NX>
struct hold_vector
{
    auto operator()(const Vector<double, NX>& x, const Vector<double, 1>&) const -> Vector<double, NX> { return x; }
};

struct constant_attitude_measurement
{
    auto operator()(const Eigen::Quaterniond&, const Vector<double, 3>&) const -> Vector<double, 1> { return Vector<double, 1>::Zero(); }

    auto jacobian(const Eigen::Quaterniond&, const Vector<double, 3>&) const -> Matrix<double, 1, 6>
    {
        Matrix<double, 1, 6> H = Matrix<double, 1, 6>::Zero();
        H(0, 0) = 1.0;
        H(0, 5) = 3.0;
        return H;
    }
};

struct hold_attitude
{
    auto operator()(const Eigen::Quaterniond& q, const Vector<double, 3>&) const -> Eigen::Quaterniond { return q; }
};

struct tilt_measurement
{
    auto operator()(const Eigen::Quaterniond& q) const -> Vector<double, 1> { return Vector<double, 1>{ctrlpp::so3::log(q)(0) + 0.25 * q.toRotationMatrix()(2, 1)}; }
};

auto divides_differently(double b, double s) -> bool
{
    return b / s != b * (1.0 / s);
}

}

TEMPLATE_TEST_CASE("one-output gain is the column-pivoting QR solve bit for bit", "[estimation][gain]", float, double)
{
    CHECK(mismatches<TestType, 1>(1) == 0);
    CHECK(mismatches<TestType, 2>(2) == 0);
    CHECK(mismatches<TestType, 3>(3) == 0);
    CHECK(mismatches<TestType, 4>(4) == 0);
    CHECK(mismatches<TestType, 6>(6) == 0);
}

TEST_CASE("kalman one-output update applies the QR solve's gain", "[estimation][gain]")
{
    SECTION("two states, several right-hand values")
    {
        ctrlpp::discrete_state_space<double, 2, 1, 1> system;
        system.A.setIdentity();
        system.B.setZero();
        system.C << 1.0, 0.0;
        system.D.setZero();
        ctrlpp::kalman_config<double, 2, 1, 1> config;
        config.P0 << 4.0, 6.0, 6.0, 10.0;
        config.x0.setZero();
        auto filter = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(system, config));
        REQUIRE(filter.update(Vector<double, 1>{1.0}).has_value());
        REQUIRE(divides_differently(6.0, 5.0));
        CHECK(same_bits(filter.state(), qr_solve(5.0, row<double, 2>{4.0, 6.0}).transpose().eval()));
    }
    SECTION("one state, one right-hand value")
    {
        ctrlpp::discrete_state_space<double, 1, 1, 1> system;
        system.A << 1.0;
        system.B << 0.0;
        system.C << 1.0;
        system.D << 0.0;
        ctrlpp::kalman_config<double, 1, 1, 1> config;
        config.R << 5.0;
        config.P0 << 9.0;
        config.x0 << 0.0;
        auto filter = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 1, 1, 1>::create(system, config));
        REQUIRE(filter.update(Vector<double, 1>{1.0}).has_value());
        REQUIRE(divides_differently(9.0, 14.0));
        CHECK(same_bits(filter.state(), qr_solve(14.0, row<double, 1>{9.0})));
    }
    SECTION("one state, zero innovation covariance keeps a zero gain")
    {
        ctrlpp::discrete_state_space<double, 1, 1, 1> system;
        system.A << 1.0;
        system.B << 0.0;
        system.C << 1.0;
        system.D << 0.0;
        ctrlpp::kalman_config<double, 1, 1, 1> config;
        config.Q.setZero();
        config.R.setZero();
        config.P0.setZero();
        config.x0 << 0.5;
        auto filter = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 1, 1, 1>::create(system, config));
        REQUIRE(qr_solve(0.0, row<double, 1>{0.0})(0) == 0.0);
        REQUIRE(filter.update(Vector<double, 1>{0.5}).has_value());
        CHECK(same_bits(filter.state(), Vector<double, 1>{0.5}));
    }
}

TEST_CASE("ekf one-output update applies the QR solve's gain", "[estimation][gain]")
{
    ctrlpp::ekf_config<double, 2, 1, 1> config;
    config.P0 << 4.0, 6.0, 6.0, 10.0;
    config.x0.setZero();
    auto filter = ctrlpp::test::constructed(ctrlpp::ekf<double, 2, 1, 1, hold_dynamics, first_state<2>>::create(hold_dynamics{}, first_state<2>{}, config));
    REQUIRE(filter.update(Vector<double, 1>{1.0}).has_value());
    REQUIRE(divides_differently(6.0, 5.0));
    CHECK(same_bits(filter.state(), qr_solve(5.0, row<double, 2>{4.0, 6.0}).transpose().eval()));
}

// Julier points with n + kappa = 4 and a Cholesky factor of small integers make
// the predicted measurement zero and S and Pxz exact.
TEST_CASE("ukf one-output update on the QR option applies the QR solve's gain", "[estimation][gain]")
{
    SECTION("two states")
    {
        ctrlpp::ukf_config<double, 2, 1, 1> config;
        config.P0 << 4.0, 6.0, 6.0, 10.0;
        config.x0.setZero();
        config.decomposition = ctrlpp::gain_decomposition::qr;
        using filter_t = ctrlpp::ukf<double, 2, 1, 1, hold_vector<2>, first_state<2>, ctrlpp::julier_sigma_points<double, 2>>;
        auto filter = ctrlpp::test::constructed(filter_t::create(hold_vector<2>{}, first_state<2>{}, config, ctrlpp::julier_options<double>{.kappa = 2.0}));
        REQUIRE(filter.update(Vector<double, 1>{1.0}).has_value());
        REQUIRE(filter.innovation()(0) == 1.0);
        REQUIRE(divides_differently(6.0, 5.0));
        CHECK(same_bits(filter.state(), qr_solve(5.0, row<double, 2>{4.0, 6.0}).transpose().eval()));
    }
    SECTION("one state")
    {
        ctrlpp::ukf_config<double, 1, 1, 1> config;
        config.R << 5.0;
        config.P0 << 9.0;
        config.x0 << 0.0;
        config.decomposition = ctrlpp::gain_decomposition::qr;
        using filter_t = ctrlpp::ukf<double, 1, 1, 1, hold_vector<1>, first_state<1>, ctrlpp::julier_sigma_points<double, 1>>;
        auto filter = ctrlpp::test::constructed(filter_t::create(hold_vector<1>{}, first_state<1>{}, config, ctrlpp::julier_options<double>{.kappa = 3.0}));
        REQUIRE(filter.update(Vector<double, 1>{1.0}).has_value());
        REQUIRE(filter.innovation()(0) == 1.0);
        REQUIRE(divides_differently(9.0, 14.0));
        CHECK(same_bits(filter.state(), qr_solve(14.0, row<double, 1>{9.0})));
    }
}

// With P0 = I and R = 0, S = H H^T and H P = H are exact; the bias reads the
// gain's last three entries.
TEST_CASE("mekf one-output update applies the QR solve's gain", "[estimation][gain]")
{
    ctrlpp::mekf_config<double, 3, 1> config;
    config.R.setZero();
    auto filter = ctrlpp::test::constructed(ctrlpp::mekf<double, 3, 1, constant_attitude_measurement>::create(constant_attitude_measurement{}, config));
    const Eigen::Quaterniond q = filter.attitude();
    const Vector<double, 3> b = filter.bias();
    REQUIRE(filter.update(Vector<double, 1>{1.0}).has_value());

    const Vector<double, 6> K = qr_solve(10.0, constant_attitude_measurement{}.jacobian(q, b)).transpose();
    const Vector<double, 3> attitude_step = K.head<3>();
    const Vector<double, 3> bias_expected = b + K.tail<3>();
    const Eigen::Quaterniond q_expected = (q * ctrlpp::so3::exp(attitude_step)).normalized();
    REQUIRE(divides_differently(3.0, 10.0));
    CHECK(same_bits(filter.bias(), bias_expected));
    CHECK(same_bits(filter.attitude().coeffs(), q_expected.coeffs()));
}

// The sigma points come from the filter's own default strategy; S and Pxz are
// accumulated in the filter's order.
TEST_CASE("manifold_ukf one-output update applies the QR solve's gain", "[estimation][gain]")
{
    ctrlpp::manifold_ukf_config<double, 1> config;
    config.R << 0.3;
    config.P0 << 0.3, 0.05, 0.0, 0.05, 0.2, 0.01, 0.0, 0.01, 0.1;
    auto filter = ctrlpp::test::constructed(ctrlpp::manifold_ukf<double, 1, hold_attitude, tilt_measurement>::create(hold_attitude{}, tilt_measurement{}, config));
    const Eigen::Quaterniond q = filter.attitude();
    const auto sigma = ctrlpp::so3_merwe_sigma_points<double>{}.generate(q, filter.covariance());
    std::array<Vector<double, 1>, ctrlpp::so3_merwe_sigma_points<double>::num_points> z_sigma;
    Vector<double, 1> z_pred = Vector<double, 1>::Zero();
    for(std::size_t i = 0; i < sigma.points.size(); ++i)
    {
        z_sigma[i] = tilt_measurement{}(sigma.points[i]);
        z_pred += sigma.Wm[i] * z_sigma[i];
    }
    Matrix<double, 1, 1> S = Matrix<double, 1, 1>::Zero();
    Matrix<double, 3, 1> Pxz = Matrix<double, 3, 1>::Zero();
    for(std::size_t i = 0; i < sigma.points.size(); ++i)
    {
        const Vector<double, 1> dz = (z_sigma[i] - z_pred).eval();
        S += sigma.Wc[i] * dz * dz.transpose();
        Eigen::Quaterniond qi = sigma.points[i];
        if(qi.dot(q) < 0.0)
            qi.coeffs() = -qi.coeffs();
        const Vector<double, 3> dx = ctrlpp::so3::log(q.conjugate() * qi);
        Pxz += sigma.Wc[i] * dx * dz.transpose();
    }
    S += config.R;

    const Vector<double, 1> z{0.4};
    REQUIRE(filter.update(z).has_value());
    const Vector<double, 3> K = qr_solve(S(0, 0), row<double, 3>{Pxz.transpose()}).transpose();
    const Vector<double, 3> delta_phi = K * (z - z_pred).eval();
    const Eigen::Quaterniond q_expected = (q * ctrlpp::so3::exp(delta_phi)).normalized();
    CHECK(same_bits(filter.innovation(), (z - z_pred).eval()));
    CHECK(same_bits(filter.attitude().coeffs(), q_expected.coeffs()));
}
