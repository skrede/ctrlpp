#include "ctrlpp/detail/covariance_ops.h"

#include "ctrlpp/estimation/detail/filter_configuration.h"

#include "ctrlpp/lie/so3.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <cmath>
#include <limits>
#include <random>
#include <cstdint>
#include <cstring>
#include <numbers>
#include <algorithm>

namespace
{

template <typename Scalar>
using matrix2 = Eigen::Matrix<Scalar, 2, 2>;

template <typename Scalar>
auto bitwise_equal(const matrix2<Scalar>& a, const matrix2<Scalar>& b) -> bool
{
    return std::memcmp(a.data(), b.data(), sizeof(Scalar) * 4) == 0;
}

template <typename Scalar>
auto pair_matrix(Scalar diagonal, Scalar upper, Scalar lower) -> matrix2<Scalar>
{
    matrix2<Scalar> M;
    M << diagonal, upper, lower, diagonal;
    return M;
}

template <typename Scalar>
auto real_quaternion(Scalar w) -> Eigen::Quaternion<Scalar>
{
    return Eigen::Quaternion<Scalar>{w, Scalar{0}, Scalar{0}, Scalar{0}};
}

template <typename Scalar>
struct unit_norm_survey
{
    std::int32_t refused;
    Scalar worst_defect;
};

template <typename Scalar, typename Producer>
auto survey(Producer produce) -> unit_norm_survey<Scalar>
{
    unit_norm_survey<Scalar> result{0, Scalar{0}};
    for(std::int32_t draw = 0; draw < 100000; ++draw)
    {
        const Eigen::Quaternion<Scalar> q = produce();
        result.worst_defect = std::max(result.worst_defect, std::abs(q.norm() - Scalar{1}));
        if(!ctrlpp::detail::checked_unit_quaternion(q).has_value())
            ++result.refused;
    }
    return result;
}

}

TEMPLATE_TEST_CASE("An exactly symmetric matrix is admitted and left bitwise untouched", "[covariance][drift]", double, float)
{
    using Scalar = TestType;
    const Scalar dm = std::numeric_limits<Scalar>::denorm_min();
    const Scalar huge = std::numeric_limits<Scalar>::max() / Scalar{2} * Scalar{1.5};
    matrix2<Scalar> M;

    SECTION("normal scale") { M << Scalar{2}, Scalar{-0.25}, Scalar{-0.25}, Scalar{1}; }

    SECTION("subnormal entries") { M << Scalar{4} * dm, dm, dm, Scalar{3} * dm; }

    SECTION("entries above half the largest finite value") { M << huge, -huge, -huge, huge; }

    SECTION("signed zeros") { M << Scalar{1}, Scalar{-0.0}, Scalar{-0.0}, Scalar{-0.0}; }

    const matrix2<Scalar> original = M;
    REQUIRE(ctrlpp::detail::symmetrize_within_drift(M));
    CHECK(bitwise_equal(M, original));
}

TEMPLATE_TEST_CASE("An asymmetry at the counted bound is drift and becomes the midpoint", "[covariance][drift]", double, float)
{
    using Scalar = TestType;
    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    const Scalar half_ulp = eps / Scalar{2};
    REQUIRE(ctrlpp::detail::rounding_bound(2, Scalar{1}, Scalar{1}) == Scalar{2} * eps);

    const Scalar upper = Scalar{0.5};
    const Scalar lower = upper + Scalar{4} * half_ulp;
    matrix2<Scalar> M = pair_matrix(Scalar{1}, upper, lower);

    REQUIRE(ctrlpp::detail::symmetrize_within_drift(M));
    const Scalar midpoint = upper + (lower - upper) / Scalar{2};
    CHECK(M(0, 1) == midpoint);
    CHECK(M(1, 0) == midpoint);
    CHECK(M(0, 0) == Scalar{1});
    CHECK(M(1, 1) == Scalar{1});
}

TEMPLATE_TEST_CASE("An asymmetry beyond the counted bound is structure and is refused unchanged", "[covariance][drift]", double, float)
{
    using Scalar = TestType;
    const Scalar half_ulp = std::numeric_limits<Scalar>::epsilon() / Scalar{2};
    const matrix2<Scalar> original = pair_matrix(Scalar{1}, Scalar{0.5}, Scalar{0.5} + Scalar{5} * half_ulp);
    matrix2<Scalar> M = original;

    CHECK_FALSE(ctrlpp::detail::symmetrize_within_drift(M));
    CHECK(bitwise_equal(M, original));

    SECTION("a counted factor of two admits the same pair")
    {
        CHECK(ctrlpp::detail::symmetrize_within_drift(M, Scalar{2}));
        CHECK(M(0, 1) == M(1, 0));
    }

    SECTION("a structurally asymmetric matrix")
    {
        matrix2<Scalar> structural = pair_matrix(Scalar{1}, Scalar{0.5}, Scalar{0.25});
        const matrix2<Scalar> kept = structural;
        CHECK_FALSE(ctrlpp::detail::symmetrize_within_drift(structural));
        CHECK(bitwise_equal(structural, kept));
    }

    SECTION("a non-finite entry")
    {
        const Scalar inf = std::numeric_limits<Scalar>::infinity();
        matrix2<Scalar> diagonal = pair_matrix(inf, Scalar{0.5}, Scalar{0.5});
        CHECK_FALSE(ctrlpp::detail::symmetrize_within_drift(diagonal));
        const matrix2<Scalar> kept = pair_matrix(Scalar{1}, inf, Scalar{0.5});
        matrix2<Scalar> off_diagonal = kept;
        CHECK_FALSE(ctrlpp::detail::symmetrize_within_drift(off_diagonal));
        CHECK(bitwise_equal(off_diagonal, kept));
    }
}

// Below max|M_ij| = denorm_min / eps the relative term underflows, so the bound
// is the absolute N * denorm_min: two subnormal units for a two-by-two.
TEMPLATE_TEST_CASE("The asymmetry bound keeps its underflow term at subnormal scale", "[covariance][drift]", double, float)
{
    using Scalar = TestType;
    const Scalar dm = std::numeric_limits<Scalar>::denorm_min();
    const Scalar diagonal = std::ldexp(dm, 20);

    matrix2<Scalar> drift = pair_matrix(diagonal, dm, Scalar{3} * dm);
    REQUIRE(ctrlpp::detail::symmetrize_within_drift(drift));
    const Scalar midpoint = dm + (Scalar{3} * dm - dm) / Scalar{2};
    CHECK(drift(0, 1) == midpoint);
    CHECK(drift(1, 0) == midpoint);

    const matrix2<Scalar> original = pair_matrix(diagonal, dm, Scalar{4} * dm);
    matrix2<Scalar> structure = original;
    CHECK_FALSE(ctrlpp::detail::symmetrize_within_drift(structure));
    CHECK(bitwise_equal(structure, original));
}

TEMPLATE_TEST_CASE("The definiteness floor is the relative bound exactly from the underflow scale up", "[covariance][floor]", double, float)
{
    using Scalar = TestType;
    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    const Scalar dm = std::numeric_limits<Scalar>::denorm_min();
    const Scalar threshold = dm / eps;

    for(const Scalar scale : {threshold, Scalar{3} * threshold, Scalar{1}, std::numeric_limits<Scalar>::max() / Scalar{4}})
        CHECK(ctrlpp::detail::rounding_bound(2, scale, Scalar{1}) == Scalar{1} * Scalar{2} * eps * scale);
    CHECK(ctrlpp::detail::rounding_bound(2, threshold / Scalar{2}, Scalar{1}) == Scalar{2} * dm);
    CHECK(ctrlpp::detail::rounding_bound(2, Scalar{0}, Scalar{1}) == Scalar{2} * dm);
}

TEMPLATE_TEST_CASE("A subnormal-scale matrix is judged against the underflow floor", "[covariance][floor]", double, float)
{
    using Scalar = TestType;
    const Scalar dm = std::numeric_limits<Scalar>::denorm_min();
    const Scalar diagonal = std::ldexp(dm, 20);

    const matrix2<Scalar> within = Eigen::Matrix<Scalar, 2, 1>{diagonal, -dm}.asDiagonal();
    CHECK(ctrlpp::detail::psd_decomposition(within).psd);

    const matrix2<Scalar> beyond = Eigen::Matrix<Scalar, 2, 1>{diagonal, Scalar{-3} * dm}.asDiagonal();
    CHECK_FALSE(ctrlpp::detail::psd_decomposition(beyond).psd);
}

TEMPLATE_TEST_CASE("A normal-scale matrix keeps the relative definiteness floor", "[covariance][floor]", double, float)
{
    using Scalar = TestType;
    const Scalar eps = std::numeric_limits<Scalar>::epsilon();

    const matrix2<Scalar> within = Eigen::Matrix<Scalar, 2, 1>{Scalar{1}, -Scalar{2} * eps}.asDiagonal();
    CHECK(ctrlpp::detail::psd_decomposition(within).psd);

    const matrix2<Scalar> beyond = Eigen::Matrix<Scalar, 2, 1>{Scalar{1}, -Scalar{3} * eps}.asDiagonal();
    CHECK_FALSE(ctrlpp::detail::psd_decomposition(beyond).psd);
}

TEMPLATE_TEST_CASE("A quaternion norm within the counted unit-norm bound is drift", "[quaternion][drift]", double, float)
{
    using Scalar = TestType;
    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    REQUIRE(ctrlpp::detail::rounding_bound(4, Scalar{1}, Scalar{1}) == Scalar{4} * eps);

    CHECK(ctrlpp::detail::checked_unit_quaternion(real_quaternion(Scalar{1})).has_value());
    CHECK(ctrlpp::detail::checked_unit_quaternion(real_quaternion(Scalar{1} + Scalar{3} * eps)).has_value());
    CHECK(ctrlpp::detail::checked_unit_quaternion(real_quaternion(Scalar{1} - Scalar{3} * eps)).has_value());
}

TEMPLATE_TEST_CASE("A quaternion norm beyond the counted unit-norm bound is refused", "[quaternion][drift]", double, float)
{
    using Scalar = TestType;
    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    const Eigen::Quaternion<Scalar> beyond = real_quaternion(Scalar{1} + Scalar{5} * eps);

    const auto refused = ctrlpp::detail::checked_unit_quaternion(beyond);
    REQUIRE_FALSE(refused.has_value());
    CHECK(refused.error() == ctrlpp::filter_error::non_unit_quaternion);

    SECTION("a counted factor of two admits the same quaternion")
    {
        CHECK(ctrlpp::detail::checked_unit_quaternion(beyond, Scalar{2}).has_value());
    }

    SECTION("a quaternion of norm two")
    {
        const auto doubled = ctrlpp::detail::checked_unit_quaternion(real_quaternion(Scalar{2}));
        REQUIRE_FALSE(doubled.has_value());
        CHECK(doubled.error() == ctrlpp::filter_error::non_unit_quaternion);
    }
}

TEMPLATE_TEST_CASE("A zero or non-finite quaternion is degenerate ahead of the unit-norm check", "[quaternion][drift]", double, float)
{
    using Scalar = TestType;
    Eigen::Quaternion<Scalar> q = real_quaternion(Scalar{0});

    SECTION("zero") {}

    SECTION("a NaN component") { q.y() = std::numeric_limits<Scalar>::quiet_NaN(); }

    SECTION("an infinite component") { q.z() = std::numeric_limits<Scalar>::infinity(); }

    const auto refused = ctrlpp::detail::checked_unit_quaternion(q, Scalar{1});
    REQUIRE_FALSE(refused.has_value());
    CHECK(refused.error() == ctrlpp::filter_error::degenerate_quaternion);
}

TEMPLATE_TEST_CASE("Every quaternion a realistic producer emits is admitted at the default factor", "[quaternion][drift]", double, float)
{
    using Scalar = TestType;
    using vector3 = Eigen::Matrix<Scalar, 3, 1>;
    using vector4 = Eigen::Matrix<Scalar, 4, 1>;
    std::mt19937 rng{1u};
    std::normal_distribution<Scalar> gaussian;
    std::uniform_real_distribution<Scalar> angle{-std::numbers::pi_v<Scalar>, std::numbers::pi_v<Scalar>};
    unit_norm_survey<Scalar> result{0, Scalar{0}};

    SECTION("normalized Gaussian 4-vectors")
    {
        result = survey<Scalar>([&] { return Eigen::Quaternion<Scalar>{vector4{gaussian(rng), gaussian(rng), gaussian(rng), gaussian(rng)}.normalized()}; });
    }

    SECTION("angle-axis conversions")
    {
        result = survey<Scalar>([&] {
            const Scalar theta = angle(rng);
            const vector3 axis = vector3{gaussian(rng), gaussian(rng), gaussian(rng)}.normalized();
            return Eigen::Quaternion<Scalar>{Eigen::AngleAxis<Scalar>{theta, axis}};
        });
    }

    SECTION("the SO(3) exponential of Gaussian rotation vectors")
    {
        result = survey<Scalar>([&] { return ctrlpp::so3::exp<Scalar>(vector3{gaussian(rng), gaussian(rng), gaussian(rng)}); });
    }

    SECTION("yaw-pitch-roll compositions")
    {
        result = survey<Scalar>([&] {
            const Scalar yaw = angle(rng);
            const Scalar pitch = angle(rng);
            const Scalar roll = angle(rng);
            return Eigen::Quaternion<Scalar>{Eigen::AngleAxis<Scalar>{yaw, vector3::UnitZ()} * Eigen::AngleAxis<Scalar>{pitch, vector3::UnitY()} * Eigen::AngleAxis<Scalar>{roll, vector3::UnitX()}};
        });
    }

    CAPTURE(result.worst_defect / std::numeric_limits<Scalar>::epsilon());
    CHECK(result.refused == 0);
}

TEST_CASE("A double unit quaternion rounded to float is admitted at the default factor", "[quaternion][drift]")
{
    std::mt19937 rng{1u};
    std::normal_distribution<double> gaussian;
    const unit_norm_survey<float> result = survey<float>([&] {
        const Eigen::Vector4d v{gaussian(rng), gaussian(rng), gaussian(rng), gaussian(rng)};
        return Eigen::Quaternionf{v.normalized().cast<float>()};
    });

    CAPTURE(result.worst_defect / std::numeric_limits<float>::epsilon());
    CHECK(result.refused == 0);
}
