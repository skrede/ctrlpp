#include "ctrlpp/nmpc.h"

#include "ctrlpp/mpc/nmpc_config.h"

#include "stub_nlp_solver.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <cstdint>
#include <utility>

namespace {

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr std::size_t NH = 5;

struct special_member_count
{
    std::int32_t defaults;
    std::int32_t copies;
    std::int32_t moves;
};

special_member_count counted{0, 0, 0};

// A by-value solver temporary shows up as a move the in-place form never makes.
struct counting_solver : ctrlpp_test::stub_nlp_solver<double>
{
    counting_solver()
    {
        ++counted.defaults;
    }

    counting_solver(const counting_solver &other)
            : stub_nlp_solver(other)
    {
        ++counted.copies;
    }

    counting_solver(counting_solver &&other) noexcept
            : stub_nlp_solver(std::move(other))
    {
        ++counted.moves;
    }

    counting_solver &operator=(const counting_solver &) = default;

    counting_solver &operator=(counting_solver &&) noexcept = default;

    ~counting_solver() = default;
};

struct double_integrator
{
    Eigen::Vector2d operator()(const Eigen::Vector2d &x, const Eigen::Matrix<double, 1, 1> &u) const
    {
        return Eigen::Vector2d{x(0) + 0.5 * x(1), x(1) + 0.5 * u(0)};
    }
};

using controller = ctrlpp::nmpc_static<double, NX, NU, NH, counting_solver, double_integrator>;

ctrlpp::nmpc_config<double, NX, NU> config()
{
    ctrlpp::nmpc_config<double, NX, NU> settings;
    settings.horizon = static_cast<int>(NH);
    return settings;
}

}

TEST_CASE("The default-solver constructor builds its solver in place", "[nmpc][static][construction]")
{
    counted = {0, 0, 0};
    controller built{double_integrator{}, config()};

    REQUIRE(counted.defaults == 1);
    REQUIRE(counted.copies == 0);
    REQUIRE(counted.moves == 0);
    REQUIRE(built.solve(Eigen::Vector2d{1.0, 0.0}).has_value());
}

TEST_CASE("The solver-taking constructor moves the caller's solver and never copies it", "[nmpc][static][construction]")
{
    counting_solver supplied;
    counted = {0, 0, 0};
    controller built{double_integrator{}, config(), std::move(supplied)};

    REQUIRE(counted.defaults == 0);
    REQUIRE(counted.copies == 0);
    REQUIRE(built.solve(Eigen::Vector2d{1.0, 0.0}).has_value());
}
