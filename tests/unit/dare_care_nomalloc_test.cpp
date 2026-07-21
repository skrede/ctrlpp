// Verify the hot paths of ctrlpp::dare and ctrlpp::care do zero heap allocation
// on fixed-size templated inputs, using the belt-and-suspenders harness from
// nomalloc_harness.h: a throwing eigen_assert that survives -DNDEBUG plus a
// global allocation counter that catches heap traffic outside Eigen's own
// bookkeeping. The harness header must stay the first include of this file.
//
// The check is surgical: we exclude the initial instantiation/warm-up and only
// guard the steady-state call itself. A negative-control case proves that both
// detection mechanisms fire, so the suite cannot silently false-pass.

#include "nomalloc_harness.h"

#include "ctrlpp/control/care.h"
#include "ctrlpp/control/dare.h"

#include "ctrlpp/detail/care_methods.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

#include <new>
#include <cstddef>
#include <utility>
#include <stdexcept>


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

// Runs the callable inside an armed no-malloc window and returns the number of
// heap allocations it performed. The count is sampled before the guard is
// released and before any test macro runs, so framework-internal allocations
// cannot pollute it.
template <typename Fn>
std::size_t guarded_allocations(Fn&& fn)
{
    ctrlpp_test::scoped_no_malloc guard;
    std::forward<Fn>(fn)();
    return guard.allocations();
}

// Warm-up-then-arm pattern shared by every steady-state case: one solve
// outside the window flushes lazy one-time instantiation, then the same solve
// must complete inside the armed window without throwing (Eigen-side trap) and
// without touching the global allocation counter.
template <typename Solve>
void require_alloc_free_steady_state(Solve&& solve)
{
    auto warmup = solve();
    REQUIRE(warmup.has_value());

    std::size_t allocations = 0;
    bool solved = false;
    allocations = guarded_allocations([&] {
        solved = solve().has_value();
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());

    REQUIRE(allocations == 0);
    REQUIRE(solved);
}

}


TEST_CASE("harness detects heap allocation",
          "[hardening][nomalloc]")
{
    SECTION("global counter fires on operator new inside the armed window")
    {
        std::size_t allocations = guarded_allocations([] {
            // Call the replaced allocation function directly: unlike a
            // new-expression, a plain function call cannot be elided.
            void* heap_block = ::operator new(sizeof(double));
            ::operator delete(heap_block);
        });

        REQUIRE(allocations > 0);
    }

    SECTION("eigen_assert sentinel fires on an Eigen allocation under -DNDEBUG")
    {
        ctrlpp_test::scoped_no_malloc guard;

        // Constructing a dynamically sized vector goes through Eigen's
        // aligned allocation check, which sets the pollable sentinel while the
        // window is armed even when the stock assert is compiled out.
        Eigen::VectorXd forced(1);
        (void)forced;

        REQUIRE(guard.eigen_violation());
    }
}

TEST_CASE("dare hot path performs zero heap allocation (NX=2, NU=1)",
          "[dare][hardening][nomalloc]")
{
    auto [A, B, Q, R] = build_dare_inputs<2, 1>();

    require_alloc_free_steady_state([&] {
        return ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    });
}

TEST_CASE("dare hot path performs zero heap allocation (NX=4, NU=2)",
          "[dare][hardening][nomalloc]")
{
    auto [A, B, Q, R] = build_dare_inputs<4, 2>();

    require_alloc_free_steady_state([&] {
        return ctrlpp::dare<double, 4, 2>(A, B, Q, R);
    });
}

TEST_CASE("dare hot path performs zero heap allocation (NX=8, NU=2)",
          "[dare][hardening][nomalloc]")
{
    auto [A, B, Q, R] = build_dare_inputs<8, 2>();

    require_alloc_free_steady_state([&] {
        return ctrlpp::dare<double, 8, 2>(A, B, Q, R);
    });
}

TEST_CASE("care hot path performs zero heap allocation (NX=2, NU=1)",
          "[care][hardening][nomalloc]")
{
    auto [A, B, Q, R] = build_care_inputs<2, 1>();

    require_alloc_free_steady_state([&] {
        return ctrlpp::care<double, 2, 1>(A, B, Q, R);
    });
}

TEST_CASE("care hot path performs zero heap allocation (NX=4, NU=2)",
          "[care][hardening][nomalloc]")
{
    auto [A, B, Q, R] = build_care_inputs<4, 2>();

    require_alloc_free_steady_state([&] {
        return ctrlpp::care<double, 4, 2>(A, B, Q, R);
    });
}

TEST_CASE("care hot path performs zero heap allocation (NX=8, NU=2)",
          "[care][hardening][nomalloc]")
{
    auto [A, B, Q, R] = build_care_inputs<8, 2>();

    require_alloc_free_steady_state([&] {
        return ctrlpp::care<double, 8, 2>(A, B, Q, R);
    });
}

TEST_CASE("care hot path performs zero heap allocation (NX=2, NU=1, sign_function)",
          "[care][hardening][nomalloc]")
{
    using ctrlpp::detail::sign_function_care_method;
    auto [A, B, Q, R] = build_care_inputs<2, 1>();

    require_alloc_free_steady_state([&] {
        return ctrlpp::care<double, 2, 1, sign_function_care_method>(A, B, Q, R);
    });
}

TEST_CASE("care hot path performs zero heap allocation (NX=4, NU=2, sign_function)",
          "[care][hardening][nomalloc]")
{
    using ctrlpp::detail::sign_function_care_method;
    auto [A, B, Q, R] = build_care_inputs<4, 2>();

    require_alloc_free_steady_state([&] {
        return ctrlpp::care<double, 4, 2, sign_function_care_method>(A, B, Q, R);
    });
}

TEST_CASE("care hot path performs zero heap allocation (NX=8, NU=2, sign_function)",
          "[care][hardening][nomalloc]")
{
    using ctrlpp::detail::sign_function_care_method;
    auto [A, B, Q, R] = build_care_inputs<8, 2>();

    require_alloc_free_steady_state([&] {
        return ctrlpp::care<double, 8, 2, sign_function_care_method>(A, B, Q, R);
    });
}

TEST_CASE("care hot path performs zero heap allocation (NX=2, NU=1, balanced_schur)",
          "[care][hardening][nomalloc]")
{
    using ctrlpp::detail::balanced_schur_care_method;
    auto [A, B, Q, R] = build_care_inputs<2, 1>();

    require_alloc_free_steady_state([&] {
        return ctrlpp::care<double, 2, 1, balanced_schur_care_method>(A, B, Q, R);
    });
}

TEST_CASE("care hot path performs zero heap allocation (NX=4, NU=2, balanced_schur)",
          "[care][hardening][nomalloc]")
{
    using ctrlpp::detail::balanced_schur_care_method;
    auto [A, B, Q, R] = build_care_inputs<4, 2>();

    require_alloc_free_steady_state([&] {
        return ctrlpp::care<double, 4, 2, balanced_schur_care_method>(A, B, Q, R);
    });
}

TEST_CASE("care hot path performs zero heap allocation (NX=8, NU=2, balanced_schur)",
          "[care][hardening][nomalloc]")
{
    using ctrlpp::detail::balanced_schur_care_method;
    auto [A, B, Q, R] = build_care_inputs<8, 2>();

    require_alloc_free_steady_state([&] {
        return ctrlpp::care<double, 8, 2, balanced_schur_care_method>(A, B, Q, R);
    });
}
