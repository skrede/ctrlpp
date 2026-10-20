#include "ctrlpp/detail/covariance_ops.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include <Eigen/Core>

#include <cmath>
#include <limits>
#include <cstring>

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
