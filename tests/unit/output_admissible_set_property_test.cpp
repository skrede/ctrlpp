#include "output_admissible_oracle.h"

#include "ctrlpp/mpc/terminal_set.h"
#include "ctrlpp/mpc/invariant_set.h"
#include "ctrlpp/mpc/output_admissible_set.h"

#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <limits>
#include <random>
#include <cstddef>
#include <cstdint>

namespace
{

namespace oracle = ctrlpp_test::oracle;

constexpr double inf = std::numeric_limits<double>::infinity();
constexpr std::int32_t explicit_samples = 4096;
constexpr std::int32_t sweep_samples = 256;

template <std::size_t N>
using point = oracle::point<N>;

template <std::size_t N>
using square = oracle::square<N>;

template <std::size_t N, std::size_t M>
struct plant
{
    square<N> A;
    ctrlpp::Matrix<double, N, M> B;
};

template <std::size_t N, std::size_t M>
struct box_bounds
{
    point<M> u_min;
    point<M> u_max;
    point<N> x_min;
    point<N> x_max;
};

template <std::size_t N, std::size_t M>
box_bounds<N, M> symmetric_bounds(double u_bound, double x_bound)
{
    return {.u_min = point<M>::Constant(-u_bound), .u_max = point<M>::Constant(u_bound), .x_min = point<N>::Constant(-x_bound), .x_max = point<N>::Constant(x_bound)};
}

template <std::size_t N, std::size_t M>
oracle::closed_loop<N> loop_of(const square<N>& A_K, const ctrlpp::Matrix<double, M, N>& K, const box_bounds<N, M>& b)
{
    return oracle::make_closed_loop<N, M>(A_K, K, b.u_min, b.u_max, b.x_min, b.x_max);
}

template <std::size_t N>
bool first_steps_span(const oracle::closed_loop<N>& loop)
{
    return oracle::positively_spans<N>(oracle::horizon_rows<N>(loop.outputs, loop.A_K, static_cast<std::int64_t>(N) - 1).H);
}

template <std::size_t N>
void require_no_violation(const ctrlpp::polytopic_set<double, N>& set, const oracle::closed_loop<N>& loop, std::mt19937& g, std::int32_t samples)
{
    REQUIRE(loop.lyapunov.decay > 0.0);
    auto const verdict = oracle::check<N>(set, loop, g, samples);
    INFO("rows " << set.H.rows() << ", vertices " << verdict.vertices << ", depth " << verdict.depth);
    CHECK_FALSE(verdict.unbounded);
    CHECK(verdict.vertex_outside == 0);
    CHECK(verdict.admissible_outside == 0);
}

template <std::size_t N>
void require_recession_refusal(const square<N>& A_K, const box_bounds<N, 1>& b)
{
    ctrlpp::Matrix<double, 1, N> const K = ctrlpp::Matrix<double, 1, N>::Zero();
    auto const result = ctrlpp::compute_output_admissible_set<double, N, 1>(A_K, K, b.u_min, b.u_max, b.x_min, b.x_max);
    CHECK_FALSE(first_steps_span<N>(loop_of<N, 1>(A_K, K, b)));
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::terminal_set_error::unbounded_terminal_set);
}

box_bounds<2, 1> unbounded_input(const point<2>& x_max)
{
    return {.u_min = point<1>::Constant(-inf), .u_max = point<1>::Constant(inf), .x_min = point<2>::Constant(-inf), .x_max = x_max};
}

template <std::size_t N, std::size_t M>
ctrlpp::Matrix<double, M, N> lqr_gain(const plant<N, M>& sys, const square<N>& P, const ctrlpp::Matrix<double, M, M>& R)
{
    ctrlpp::Matrix<double, M, M> const RpBtPB = R + sys.B.transpose() * P * sys.B;
    return -(RpBtPB.ldlt().solve(sys.B.transpose() * P * sys.A));
}

plant<2, 1> double_integrator()
{
    constexpr double dt = 0.1;
    return {.A = (square<2>() << 1.0, dt, 0.0, 1.0).finished(), .B = (ctrlpp::Matrix<double, 2, 1>() << 0.5 * dt * dt, dt).finished()};
}

struct sweep_counts
{
    std::int32_t returned = 0;
    std::int32_t ingredients_refused = 0;
    std::int32_t unbounded = 0;
    std::int32_t unbounded_although_spanning = 0;
    std::int32_t not_converged = 0;
    std::int32_t other_refusal = 0;
};

void tally_refusal(sweep_counts& counts, ctrlpp::terminal_set_error error)
{
    if(error == ctrlpp::terminal_set_error::unbounded_terminal_set)
        ++counts.unbounded;
    else if(error == ctrlpp::terminal_set_error::not_converged)
        ++counts.not_converged;
    else
        ++counts.other_refusal;
}

// One plant through terminal_ingredients with Q = I and R = 1. A returned set must pass the oracle. A refusal is a
// defect when the bounds are two-sided, since every hypothesis then holds; with one-sided bounds the refusal must be
// the recession one whenever the rows of the first N steps fail to positively span.
template <std::size_t N>
void run_case(const plant<N, 1>& sys, const box_bounds<N, 1>& b, bool two_sided, std::mt19937& sampler, sweep_counts& counts)
{
    ctrlpp::Matrix<double, 1, 1> const R = ctrlpp::Matrix<double, 1, 1>::Identity();
    auto const ingredients = ctrlpp::terminal_ingredients<double, N, 1>(sys.A, sys.B, square<N>::Identity(), R, b.u_min, b.u_max, b.x_min, b.x_max);
    if(!ingredients)
    {
        ++counts.ingredients_refused;
        return;
    }
    auto const K = lqr_gain<N, 1>(sys, ingredients->Qf, R);
    auto const loop = loop_of<N, 1>(sys.A + sys.B * K, K, b);
    INFO("A =\n" << sys.A << "\nB = " << sys.B.transpose() << "\nx_min = " << b.x_min.transpose() << "\nx_max = " << b.x_max.transpose());
    if(!ingredients->invariant_set)
    {
        bool const spans = first_steps_span<N>(loop);
        tally_refusal(counts, ingredients->invariant_set.error());
        counts.unbounded_although_spanning += spans && ingredients->invariant_set.error() == ctrlpp::terminal_set_error::unbounded_terminal_set ? 1 : 0;
        CHECK_FALSE(two_sided);
        CHECK((spans || ingredients->invariant_set.error() == ctrlpp::terminal_set_error::unbounded_terminal_set));
        return;
    }
    ++counts.returned;
    CHECK(first_steps_span<N>(loop));
    require_no_violation<N>(*ingredients->invariant_set, loop, sampler, sweep_samples);
}

template <std::size_t N>
plant<N, 1> generic_plant(std::mt19937& g)
{
    plant<N, 1> sys;
    for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(N); ++i)
    {
        for(Eigen::Index j = 0; j < static_cast<Eigen::Index>(N); ++j)
            sys.A(i, j) = 2.0 * oracle::unit_draw(g) - 1.0;
        sys.B(i) = 2.0 * oracle::unit_draw(g) - 1.0;
    }
    return sys;
}

template <std::size_t N>
plant<N, 1> half_integer_plant(std::mt19937& g)
{
    plant<N, 1> sys;
    for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(N); ++i)
    {
        for(Eigen::Index j = 0; j < static_cast<Eigen::Index>(N); ++j)
            sys.A(i, j) = oracle::half_integer_draw(g);
        sys.B(i) = oracle::half_integer_draw(g);
    }
    return sys;
}

template <std::size_t N>
box_bounds<N, 1> asymmetric_box(std::mt19937& g)
{
    auto b = symmetric_bounds<N, 1>(1.0, 1.0);
    for(Eigen::Index j = 0; j < static_cast<Eigen::Index>(N); ++j)
    {
        b.x_min(j) = -0.5 - oracle::unit_draw(g);
        b.x_max(j) = 0.5 + oracle::unit_draw(g);
    }
    return b;
}

template <std::size_t N>
box_bounds<N, 1> one_sided_box(std::mt19937& g, bool inputs_too)
{
    auto b = asymmetric_box<N>(g);
    for(Eigen::Index j = 0; j < static_cast<Eigen::Index>(N); ++j)
    {
        if(g() % 2U == 0U)
            b.x_min(j) = -inf;
        if(g() % 2U == 0U)
            b.x_max(j) = inf;
    }
    if(inputs_too && g() % 2U == 0U)
        b.u_min(0) = -inf;
    if(inputs_too && g() % 2U == 0U)
        b.u_max(0) = inf;
    return b;
}

void report(const char* name, const sweep_counts& counts)
{
    WARN(name << ": returned " << counts.returned << ", unbounded_terminal_set " << counts.unbounded << " (" << counts.unbounded_although_spanning << " although the first steps span), not_converged "
              << counts.not_converged << ", other refusals " << counts.other_refusal << ", terminal_ingredients refusals " << counts.ingredients_refused);
}

template <std::size_t N>
void one_sided_sweep(std::uint32_t seed, std::int32_t cases, bool inputs_too)
{
    std::mt19937 g(seed);
    std::mt19937 sampler(g());
    sweep_counts counts;
    for(std::int32_t c = 0; c < cases; ++c)
    {
        auto const sys = generic_plant<N>(g);
        run_case<N>(sys, one_sided_box<N>(g, inputs_too), false, sampler, counts);
    }
    report(inputs_too ? "one-sided states and inputs" : "one-sided states", counts);
}

template <std::size_t N>
void two_sided_sweep(const char* name, std::uint32_t seed, std::int32_t cases, bool half_integer)
{
    std::mt19937 g(seed);
    std::mt19937 sampler(g());
    sweep_counts counts;
    for(std::int32_t c = 0; c < cases; ++c)
    {
        auto const sys = half_integer ? half_integer_plant<N>(g) : generic_plant<N>(g);
        auto const b = half_integer ? symmetric_bounds<N, 1>(1.0, 1.0) : asymmetric_box<N>(g);
        run_case<N>(sys, b, true, sampler, counts);
    }
    report(name, counts);
    CHECK(counts.unbounded + counts.not_converged + counts.other_refusal == 0);
}

}

TEST_CASE("A single one-sided state bound leaves a recession direction and is refused", "[output_admissible][property][one_sided]")
{
    square<2> const A_K = (square<2>() << 0.5, 1.0, 0.0, 0.5).finished();
    require_recession_refusal<2>(A_K, unbounded_input(point<2>(1.0, inf)));
}

TEST_CASE("Upper bounds alone on a coupled loop leave a recession direction and are refused", "[output_admissible][property][one_sided]")
{
    square<2> const A_K = (square<2>() << 0.5, 1.0, 0.0, 0.5).finished();
    require_recession_refusal<2>(A_K, unbounded_input(point<2>(1.0, 1.0)));
}

TEST_CASE("Upper bounds alone on a rotating loop are refused although later steps would close the set", "[output_admissible][property][one_sided]")
{
    square<2> const A_K = (square<2>() << 0.9, 0.3, -0.3, 0.9).finished();
    require_recession_refusal<2>(A_K, unbounded_input(point<2>(1.0, 1.0)));
}

TEST_CASE("A one-sided state bound with two-sided input bounds yields a bounded admissible set", "[output_admissible][property][one_sided]")
{
    auto const sys = double_integrator();
    box_bounds<2, 1> const b{.u_min = point<1>::Constant(-1.0), .u_max = point<1>::Constant(1.0), .x_min = point<2>::Constant(-inf), .x_max = point<2>(10.0, inf)};
    ctrlpp::Matrix<double, 1, 1> const R = ctrlpp::Matrix<double, 1, 1>::Identity();
    auto const ingredients = ctrlpp::terminal_ingredients<double, 2, 1>(sys.A, sys.B, square<2>::Identity(), R, b.u_min, b.u_max, b.x_min, b.x_max);
    REQUIRE(ingredients.has_value());
    REQUIRE(ingredients->invariant_set.has_value());
    auto const K = lqr_gain<2, 1>(sys, ingredients->Qf, R);
    auto const loop = loop_of<2, 1>(sys.A + sys.B * K, K, b);
    CHECK(first_steps_span<2>(loop));
    std::mt19937 g(1U);
    require_no_violation<2>(*ingredients->invariant_set, loop, g, explicit_samples);
}

TEST_CASE("One-sided state bounds at two states", "[output_admissible][property][one_sided][sweep]")
{
    one_sided_sweep<2>(21U, 1000, false);
}

TEST_CASE("One-sided state bounds at three states", "[output_admissible][property][one_sided][sweep]")
{
    one_sided_sweep<3>(31U, 500, false);
}

TEST_CASE("One-sided state bounds at four states", "[output_admissible][property][one_sided][sweep]")
{
    one_sided_sweep<4>(41U, 200, false);
}

TEST_CASE("One-sided state and input bounds at two states", "[output_admissible][property][one_sided][sweep]")
{
    one_sided_sweep<2>(23U, 1000, true);
}

TEST_CASE("One-sided state and input bounds at three states", "[output_admissible][property][one_sided][sweep]")
{
    one_sided_sweep<3>(34U, 500, true);
}

TEST_CASE("One-sided state and input bounds at four states", "[output_admissible][property][one_sided][sweep]")
{
    one_sided_sweep<4>(43U, 200, true);
}

TEST_CASE("A rotation whose degenerate corners repeat vertices yields a bounded admissible set", "[output_admissible][property][degenerate]")
{
    plant<2, 1> const sys{.A = (square<2>() << 0.0, -1.0, 1.0, 0.0).finished(), .B = (ctrlpp::Matrix<double, 2, 1>() << -1.0, 1.0).finished()};
    auto const b = symmetric_bounds<2, 1>(1.0, 1.0);
    ctrlpp::Matrix<double, 1, 1> const R = ctrlpp::Matrix<double, 1, 1>::Identity();
    auto const ingredients = ctrlpp::terminal_ingredients<double, 2, 1>(sys.A, sys.B, square<2>::Identity(), R, b.u_min, b.u_max, b.x_min, b.x_max);
    REQUIRE(ingredients.has_value());
    REQUIRE(ingredients->invariant_set.has_value());
    auto const K = lqr_gain<2, 1>(sys, ingredients->Qf, R);
    std::mt19937 g(2U);
    require_no_violation<2>(*ingredients->invariant_set, loop_of<2, 1>(sys.A + sys.B * K, K, b), g, explicit_samples);
}

TEST_CASE("An input bound coinciding with a state bound yields a bounded admissible set", "[output_admissible][property][degenerate]")
{
    square<2> const A_K = (square<2>() << 0.5, 0.4, -0.3, 0.6).finished();
    auto const gain = GENERATE(0.5, 0.3);
    ctrlpp::Matrix<double, 1, 2> const K(-gain, 0.0);
    auto const b = symmetric_bounds<2, 1>(gain, 1.0);
    auto const result = ctrlpp::compute_output_admissible_set<double, 2, 1>(A_K, K, b.u_min, b.u_max, b.x_min, b.x_max);
    INFO("gain " << gain);
    REQUIRE(result.has_value());
    std::mt19937 g(3U);
    require_no_violation<2>(*result, loop_of<2, 1>(A_K, K, b), g, explicit_samples);
}

TEST_CASE("Half-integer plants at two states", "[output_admissible][property][sweep]")
{
    two_sided_sweep<2>("two states, half-integer", 22U, 2000, true);
}

TEST_CASE("Generic plants with asymmetric boxes at three states", "[output_admissible][property][sweep]")
{
    two_sided_sweep<3>("three states, generic", 32U, 1000, false);
}

TEST_CASE("Half-integer plants at three states", "[output_admissible][property][sweep]")
{
    two_sided_sweep<3>("three states, half-integer", 33U, 1000, true);
}

TEST_CASE("Generic plants with asymmetric boxes at four states", "[output_admissible][property][sweep]")
{
    two_sided_sweep<4>("four states, generic", 42U, 100, false);
}
