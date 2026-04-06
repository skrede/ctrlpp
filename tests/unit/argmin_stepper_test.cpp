#include "ctrlpp/mpc/argmin_solver.h"
#include "ctrlpp/mpc/argmin_policies.h"
#include "ctrlpp/mpc/nlp_solver.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>

#include <cstddef>

namespace
{

auto make_rosenbrock() -> ctrlpp::nlp_problem<double>
{
    ctrlpp::nlp_problem<double> prob;
    prob.n_vars = 2;
    prob.n_constraints = 0;
    prob.cost = [](std::span<const double> x) {
        return (1.0 - x[0]) * (1.0 - x[0]) + 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]);
    };
    prob.gradient = [](std::span<const double> x, std::span<double> g) {
        g[0] = -2.0 * (1.0 - x[0]) - 400.0 * x[0] * (x[1] - x[0] * x[0]);
        g[1] = 200.0 * (x[1] - x[0] * x[0]);
    };
    prob.x_lower = Eigen::Vector2d::Constant(-10.0);
    prob.x_upper = Eigen::Vector2d::Constant(10.0);
    prob.c_lower = Eigen::VectorXd{};
    prob.c_upper = Eigen::VectorXd{};
    return prob;
}

}

TEST_CASE("argmin step concept satisfaction", "[argmin][stepper]")
{
    static_assert(ctrlpp::nlp_stepper<ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp>>,
                  "argmin_solver must satisfy nlp_stepper concept");
}

TEST_CASE("argmin step budget exhaustion", "[argmin][stepper]")
{
    auto prob = make_rosenbrock();

    ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp> solver;
    solver.setup(prob);

    ctrlpp::nlp_update<double> update;
    update.x0 = Eigen::Vector2d{-2.0, 2.0};

    auto result = solver.step(update, 3);

    CHECK(result.status == ctrlpp::solve_status::max_iterations);
    CHECK(result.iterations <= 3);
}

TEST_CASE("argmin step then solve converges", "[argmin][stepper]")
{
    auto prob = make_rosenbrock();

    ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp> solver;
    solver.setup(prob);

    ctrlpp::nlp_update<double> update;
    update.x0 = Eigen::Vector2d{-2.0, 2.0};

    // Partial solve with small budget
    auto partial = solver.step(update, 3);
    CHECK(partial.iterations <= 3);

    // Continue from partial result
    ctrlpp::nlp_update<double> continuation;
    continuation.x0 = partial.x;

    auto final_result = solver.solve(continuation);

    CHECK(final_result.status == ctrlpp::solve_status::optimal);
    CHECK(final_result.x(0) == Catch::Approx(1.0).margin(1e-3));
    CHECK(final_result.x(1) == Catch::Approx(1.0).margin(1e-3));
}
