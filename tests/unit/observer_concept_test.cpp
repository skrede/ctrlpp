#include "ctrlpp/estimation/observer_policy.h"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <type_traits>

// --- Mock observers for concept testing ---

struct MockMinimalObserver
{
    using observer_tag = void;
    using state_vector_t = std::array<double, 2>;
    using input_vector_t = std::array<double, 1>;
    using output_vector_t = std::array<double, 1>;

    void predict(const input_vector_t&) {}
    void update(const output_vector_t&) {}
    auto state() const -> const state_vector_t& { return state_; }

private:
    state_vector_t state_{};
};

struct MockCovarianceObserver
{
    using observer_tag = void;
    using state_vector_t = std::array<double, 2>;
    using input_vector_t = std::array<double, 1>;
    using output_vector_t = std::array<double, 1>;
    using covariance_t = std::array<std::array<double, 2>, 2>;
    using innovation_t = std::array<double, 1>;

    void predict(const input_vector_t&) {}
    void update(const output_vector_t&) {}
    auto state() const -> const state_vector_t& { return state_; }
    auto covariance() const -> const covariance_t& { return cov_; }
    auto innovation() const -> const innovation_t& { return innov_; }

private:
    state_vector_t state_{};
    covariance_t cov_{};
    innovation_t innov_{};
};

// --- Compile-time concept checks ---

static_assert(ctrlpp::ObserverPolicy<ctrlpp::null_observer>);
static_assert(ctrlpp::ObserverPolicy<MockMinimalObserver>);
static_assert(ctrlpp::ObserverPolicy<MockCovarianceObserver>);

static_assert(ctrlpp::CovarianceObserver<MockCovarianceObserver>);
static_assert(!ctrlpp::CovarianceObserver<MockMinimalObserver>);
static_assert(!ctrlpp::CovarianceObserver<ctrlpp::null_observer>);

static_assert(!ctrlpp::ObserverPolicy<int>);
static_assert(!ctrlpp::ObserverPolicy<double>);

// --- Runtime tests ---

TEST_CASE("null_observer predict/update/state do not crash", "[observer]")
{
    ctrlpp::null_observer obs;
    std::monostate m;

    obs.predict(m);
    obs.update(m);

    [[maybe_unused]] auto& s = obs.state();
    REQUIRE(std::is_same_v<std::remove_cvref_t<decltype(s)>, std::monostate>);
}

TEST_CASE("MockMinimalObserver satisfies ObserverPolicy", "[observer]")
{
    MockMinimalObserver obs;
    MockMinimalObserver::input_vector_t u{1.0};
    MockMinimalObserver::output_vector_t z{2.0};

    obs.predict(u);
    obs.update(z);

    auto& s = obs.state();
    REQUIRE(s[0] == 0.0);
    REQUIRE(s[1] == 0.0);
}

TEST_CASE("MockCovarianceObserver satisfies CovarianceObserver", "[observer]")
{
    MockCovarianceObserver obs;
    MockCovarianceObserver::input_vector_t u{1.0};
    MockCovarianceObserver::output_vector_t z{2.0};

    obs.predict(u);
    obs.update(z);

    auto& s = obs.state();
    auto& cov = obs.covariance();
    auto& innov = obs.innovation();

    REQUIRE(s[0] == 0.0);
    REQUIRE(cov[0][0] == 0.0);
    REQUIRE(innov[0] == 0.0);
}
