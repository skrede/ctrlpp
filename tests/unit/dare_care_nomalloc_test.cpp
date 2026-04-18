// Verify the hot paths of ctrlpp::dare and ctrlpp::care do zero heap allocation
// on fixed-size templated inputs. Uses Eigen's EIGEN_RUNTIME_NO_MALLOC + assert
// contract: any heap alloc inside a `set_is_malloc_allowed(false)` window asserts.
//
// The check is surgical: we exclude the initial instantiation/warm-up and only
// guard the steady-state call itself.

#define EIGEN_RUNTIME_NO_MALLOC

#include "ctrlpp/control/dare.h"
#include "ctrlpp/control/care.h"
#include "ctrlpp/detail/care_methods.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>


namespace
{

template <int NX, int NU>
auto build_dare_inputs()
{
    Eigen::Matrix<double, NX, NX> A = Eigen::Matrix<double, NX, NX>::Identity();
    for(int i = 0; i + 1 < NX; ++i)
        A(i, i + 1) = 0.05;

    Eigen::Matrix<double, NX, NU> B = Eigen::Matrix<double, NX, NU>::Zero();
    for(int j = 0; j < NU; ++j)
        B(std::min((j + 1) * (NX / NU), NX) - 1, j) = 0.05;

    Eigen::Matrix<double, NX, NX> Q = Eigen::Matrix<double, NX, NX>::Identity();
    Eigen::Matrix<double, NU, NU> R = 0.1 * Eigen::Matrix<double, NU, NU>::Identity();

    return std::tuple{A, B, Q, R};
}

template <int NX, int NU>
auto build_care_inputs()
{
    Eigen::Matrix<double, NX, NX> A = Eigen::Matrix<double, NX, NX>::Zero();
    for(int i = 0; i < NX; ++i)
        A(i, i) = -0.5;
    for(int i = 0; i + 1 < NX; ++i)
        A(i, i + 1) = 1.0;

    Eigen::Matrix<double, NX, NU> B = Eigen::Matrix<double, NX, NU>::Zero();
    for(int j = 0; j < NU; ++j)
        B(std::min((j + 1) * (NX / NU), NX) - 1, j) = 1.0;

    Eigen::Matrix<double, NX, NX> Q = Eigen::Matrix<double, NX, NX>::Identity();
    Eigen::Matrix<double, NU, NU> R = 0.1 * Eigen::Matrix<double, NU, NU>::Identity();

    return std::tuple{A, B, Q, R};
}

}


TEST_CASE("dare hot path performs zero heap allocation (NX=2, NU=1)",
          "[dare][hardening][nomalloc]")
{
    auto [A, B, Q, R] = build_dare_inputs<2, 1>();

    // Warm-up call outside the no-malloc window (instantiates any lazy members).
    auto warmup = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(warmup.has_value());

    Eigen::internal::set_is_malloc_allowed(false);
    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    Eigen::internal::set_is_malloc_allowed(true);

    REQUIRE(result.has_value());
}

TEST_CASE("dare hot path performs zero heap allocation (NX=4, NU=2)",
          "[dare][hardening][nomalloc]")
{
    auto [A, B, Q, R] = build_dare_inputs<4, 2>();

    auto warmup = ctrlpp::dare<double, 4, 2>(A, B, Q, R);
    REQUIRE(warmup.has_value());

    Eigen::internal::set_is_malloc_allowed(false);
    auto result = ctrlpp::dare<double, 4, 2>(A, B, Q, R);
    Eigen::internal::set_is_malloc_allowed(true);

    REQUIRE(result.has_value());
}

TEST_CASE("dare hot path performs zero heap allocation (NX=8, NU=2)",
          "[dare][hardening][nomalloc]")
{
    auto [A, B, Q, R] = build_dare_inputs<8, 2>();

    auto warmup = ctrlpp::dare<double, 8, 2>(A, B, Q, R);
    REQUIRE(warmup.has_value());

    Eigen::internal::set_is_malloc_allowed(false);
    auto result = ctrlpp::dare<double, 8, 2>(A, B, Q, R);
    Eigen::internal::set_is_malloc_allowed(true);

    REQUIRE(result.has_value());
}

TEST_CASE("care hot path performs zero heap allocation (NX=2, NU=1)",
          "[care][hardening][nomalloc]")
{
    auto [A, B, Q, R] = build_care_inputs<2, 1>();

    auto warmup = ctrlpp::care<double, 2, 1>(A, B, Q, R);
    REQUIRE(warmup.has_value());

    Eigen::internal::set_is_malloc_allowed(false);
    auto result = ctrlpp::care<double, 2, 1>(A, B, Q, R);
    Eigen::internal::set_is_malloc_allowed(true);

    REQUIRE(result.has_value());
}

TEST_CASE("care hot path performs zero heap allocation (NX=4, NU=2)",
          "[care][hardening][nomalloc]")
{
    auto [A, B, Q, R] = build_care_inputs<4, 2>();

    auto warmup = ctrlpp::care<double, 4, 2>(A, B, Q, R);
    REQUIRE(warmup.has_value());

    Eigen::internal::set_is_malloc_allowed(false);
    auto result = ctrlpp::care<double, 4, 2>(A, B, Q, R);
    Eigen::internal::set_is_malloc_allowed(true);

    REQUIRE(result.has_value());
}

TEST_CASE("care hot path performs zero heap allocation (NX=8, NU=2)",
          "[care][hardening][nomalloc]")
{
    auto [A, B, Q, R] = build_care_inputs<8, 2>();

    auto warmup = ctrlpp::care<double, 8, 2>(A, B, Q, R);
    REQUIRE(warmup.has_value());

    Eigen::internal::set_is_malloc_allowed(false);
    auto result = ctrlpp::care<double, 8, 2>(A, B, Q, R);
    Eigen::internal::set_is_malloc_allowed(true);

    REQUIRE(result.has_value());
}

TEST_CASE("care hot path performs zero heap allocation (NX=2, NU=1, sign_function)",
          "[care][hardening][nomalloc]")
{
    using ctrlpp::detail::sign_function_care_method;
    auto [A, B, Q, R] = build_care_inputs<2, 1>();

    auto warmup = ctrlpp::care<double, 2, 1, sign_function_care_method>(A, B, Q, R);
    REQUIRE(warmup.has_value());

    Eigen::internal::set_is_malloc_allowed(false);
    auto result = ctrlpp::care<double, 2, 1, sign_function_care_method>(A, B, Q, R);
    Eigen::internal::set_is_malloc_allowed(true);

    REQUIRE(result.has_value());
}

TEST_CASE("care hot path performs zero heap allocation (NX=4, NU=2, sign_function)",
          "[care][hardening][nomalloc]")
{
    using ctrlpp::detail::sign_function_care_method;
    auto [A, B, Q, R] = build_care_inputs<4, 2>();

    auto warmup = ctrlpp::care<double, 4, 2, sign_function_care_method>(A, B, Q, R);
    REQUIRE(warmup.has_value());

    Eigen::internal::set_is_malloc_allowed(false);
    auto result = ctrlpp::care<double, 4, 2, sign_function_care_method>(A, B, Q, R);
    Eigen::internal::set_is_malloc_allowed(true);

    REQUIRE(result.has_value());
}

TEST_CASE("care hot path performs zero heap allocation (NX=8, NU=2, sign_function)",
          "[care][hardening][nomalloc]")
{
    using ctrlpp::detail::sign_function_care_method;
    auto [A, B, Q, R] = build_care_inputs<8, 2>();

    auto warmup = ctrlpp::care<double, 8, 2, sign_function_care_method>(A, B, Q, R);
    REQUIRE(warmup.has_value());

    Eigen::internal::set_is_malloc_allowed(false);
    auto result = ctrlpp::care<double, 8, 2, sign_function_care_method>(A, B, Q, R);
    Eigen::internal::set_is_malloc_allowed(true);

    REQUIRE(result.has_value());
}

TEST_CASE("care hot path performs zero heap allocation (NX=2, NU=1, balanced_schur)",
          "[care][hardening][nomalloc]")
{
    using ctrlpp::detail::balanced_schur_care_method;
    auto [A, B, Q, R] = build_care_inputs<2, 1>();

    auto warmup = ctrlpp::care<double, 2, 1, balanced_schur_care_method>(A, B, Q, R);
    REQUIRE(warmup.has_value());

    Eigen::internal::set_is_malloc_allowed(false);
    auto result = ctrlpp::care<double, 2, 1, balanced_schur_care_method>(A, B, Q, R);
    Eigen::internal::set_is_malloc_allowed(true);

    REQUIRE(result.has_value());
}

TEST_CASE("care hot path performs zero heap allocation (NX=4, NU=2, balanced_schur)",
          "[care][hardening][nomalloc]")
{
    using ctrlpp::detail::balanced_schur_care_method;
    auto [A, B, Q, R] = build_care_inputs<4, 2>();

    auto warmup = ctrlpp::care<double, 4, 2, balanced_schur_care_method>(A, B, Q, R);
    REQUIRE(warmup.has_value());

    Eigen::internal::set_is_malloc_allowed(false);
    auto result = ctrlpp::care<double, 4, 2, balanced_schur_care_method>(A, B, Q, R);
    Eigen::internal::set_is_malloc_allowed(true);

    REQUIRE(result.has_value());
}

TEST_CASE("care hot path performs zero heap allocation (NX=8, NU=2, balanced_schur)",
          "[care][hardening][nomalloc]")
{
    using ctrlpp::detail::balanced_schur_care_method;
    auto [A, B, Q, R] = build_care_inputs<8, 2>();

    auto warmup = ctrlpp::care<double, 8, 2, balanced_schur_care_method>(A, B, Q, R);
    REQUIRE(warmup.has_value());

    Eigen::internal::set_is_malloc_allowed(false);
    auto result = ctrlpp::care<double, 8, 2, balanced_schur_care_method>(A, B, Q, R);
    Eigen::internal::set_is_malloc_allowed(true);

    REQUIRE(result.has_value());
}
