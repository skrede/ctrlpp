#include "ctrlpp/estimation/resampling/multinomial_resampling.h"
#include "ctrlpp/estimation/resampling/resampling_strategy.h"
#include "ctrlpp/estimation/resampling/systematic_resampling.h"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstddef>
#include <random>

using namespace ctrlpp;

// ---------------------------------------------------------------------------
// Concept satisfaction: resampling strategies
// ---------------------------------------------------------------------------
static_assert(resampling_strategy<systematic_resampling, std::mt19937_64, 10>);
static_assert(resampling_strategy<multinomial_resampling, std::mt19937_64, 10>);

// ---------------------------------------------------------------------------
// Systematic resampling tests
// ---------------------------------------------------------------------------

TEST_CASE("systematic resampling uniform weights produce identity indices")
{
    constexpr std::size_t NP = 8;
    std::array<double, NP> weights;
    weights.fill(1.0 / static_cast<double>(NP));

    std::array<std::size_t, NP> indices{};
    std::mt19937_64 rng(42);

    systematic_resampling resampler;
    resampler.resample(weights, indices, rng);

    for(std::size_t i = 0; i < NP; ++i)
    {
        CHECK(indices[i] == i);
    }
}

TEST_CASE("systematic resampling degenerate weights concentrate on single index")
{
    constexpr std::size_t NP = 10;
    std::array<double, NP> weights{};
    weights.fill(0.0);
    weights[3] = 1.0;

    std::array<std::size_t, NP> indices{};
    std::mt19937_64 rng(123);

    systematic_resampling resampler;
    resampler.resample(weights, indices, rng);

    for(std::size_t i = 0; i < NP; ++i)
    {
        CHECK(indices[i] == 3);
    }
}

TEST_CASE("systematic resampling all indices valid")
{
    constexpr std::size_t NP = 100;
    std::array<double, NP> weights;
    weights.fill(1.0 / static_cast<double>(NP));

    std::array<std::size_t, NP> indices{};
    std::mt19937_64 rng(7);

    systematic_resampling resampler;
    resampler.resample(weights, indices, rng);

    for(std::size_t i = 0; i < NP; ++i)
    {
        CHECK(indices[i] < NP);
    }
}

TEST_CASE("systematic resampling deterministic with fixed seed")
{
    constexpr std::size_t NP = 20;
    std::array<double, NP> weights;
    for(std::size_t i = 0; i < NP; ++i)
    {
        weights[i] = static_cast<double>(i + 1);
    }
    double sum = 0.0;
    for(auto w : weights)
        sum += w;
    for(auto& w : weights)
        w /= sum;

    std::array<std::size_t, NP> indices1{};
    std::array<std::size_t, NP> indices2{};

    {
        std::mt19937_64 rng(999);
        systematic_resampling resampler;
        resampler.resample(weights, indices1, rng);
    }
    {
        std::mt19937_64 rng(999);
        systematic_resampling resampler;
        resampler.resample(weights, indices2, rng);
    }

    for(std::size_t i = 0; i < NP; ++i)
    {
        CHECK(indices1[i] == indices2[i]);
    }
}

// ---------------------------------------------------------------------------
// Multinomial resampling tests
// ---------------------------------------------------------------------------

TEST_CASE("multinomial resampling degenerate weights concentrate on single index")
{
    constexpr std::size_t NP = 10;
    std::array<double, NP> weights{};
    weights.fill(0.0);
    weights[5] = 1.0;

    std::array<std::size_t, NP> indices{};
    std::mt19937_64 rng(456);

    multinomial_resampling resampler;
    resampler.resample(weights, indices, rng);

    for(std::size_t i = 0; i < NP; ++i)
    {
        CHECK(indices[i] == 5);
    }
}

TEST_CASE("multinomial resampling uniform weights roughly uniform distribution")
{
    constexpr std::size_t NP = 1000;
    std::array<double, NP> weights;
    weights.fill(1.0 / static_cast<double>(NP));

    std::array<std::size_t, NP> indices{};
    std::mt19937_64 rng(789);

    multinomial_resampling resampler;
    resampler.resample(weights, indices, rng);

    for(std::size_t i = 0; i < NP; ++i)
    {
        CHECK(indices[i] < NP);
    }

    std::array<int, NP> counts{};
    for(auto idx : indices)
    {
        ++counts[idx];
    }
    std::size_t unique = 0;
    for(auto c : counts)
    {
        if(c > 0)
            ++unique;
    }
    CHECK(unique > NP / 2);
}
