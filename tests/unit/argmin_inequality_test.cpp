#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/argmin_solver.h"
#include "ctrlpp/mpc/argmin_policies.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>

#include <span>
#include <cmath>
#include <limits>

namespace
{

constexpr double inf = std::numeric_limits<double>::infinity();

using solver_type = ctrlpp::argmin_solver<double, ctrlpp::argmin_nw_sqp>;

enum class jacobian_path
{
    analytic,
    finite_difference
};

// Minimizes 0.5 |x - (1, 1)|^2 subject to lower <= x0 + x1 <= upper, posed as one general row.
auto make_problem(double lower, double upper, jacobian_path path) -> ctrlpp::nlp_problem<double>
{
    ctrlpp::nlp_problem<double> problem{};
    problem.n_vars        = 2;
    problem.n_constraints = 1;
    problem.cost          = [](std::span<const double> x) { return 0.5 * ((x[0] - 1.0) * (x[0] - 1.0) + (x[1] - 1.0) * (x[1] - 1.0)); };
    problem.gradient      = [](std::span<const double> x, std::span<double> g)
    {
        g[0] = x[0] - 1.0;
        g[1] = x[1] - 1.0;
    };
    problem.constraints = [](std::span<const double> x, std::span<double> c) { c[0] = x[0] + x[1]; };
    if(path == jacobian_path::analytic)
    {
        problem.constraint_jacobian = [](std::span<const double>, std::span<double> jac)
        {
            jac[0] = 1.0;
            jac[1] = 1.0;
        };
    }
    problem.x_lower = Eigen::Vector2d::Constant(-10.0);
    problem.x_upper = Eigen::Vector2d::Constant(10.0);
    problem.c_lower = Eigen::VectorXd::Constant(1, lower);
    problem.c_upper = Eigen::VectorXd::Constant(1, upper);
    return problem;
}

// On the active set, x - x* = K^-1 r for the KKT matrix K of the identity Hessian and the active rows, and kkt_tol
// bounds every entry of the residual r, so |x - x*| <= sqrt(dim r) * kkt_tol / sigma_min(K). The finite-difference
// Jacobian of a linear row is exact to rounding, orders of magnitude below kkt_tol.
auto optimum_tolerance(Eigen::Index active_rows) -> double
{
    const Eigen::Index n = 2;
    Eigen::MatrixXd kkt  = Eigen::MatrixXd::Zero(n + active_rows, n + active_rows);
    kkt.topLeftCorner(n, n).setIdentity();
    kkt.topRightCorner(n, active_rows).setOnes();
    kkt.bottomLeftCorner(active_rows, n).setOnes();
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(kkt);
    return std::sqrt(static_cast<double>(n + active_rows)) * ctrlpp::argmin_settings<double>{}.kkt_tol / svd.singularValues().minCoeff();
}

void check_path(const ctrlpp::nlp_problem<double>& problem, double optimum_sum, Eigen::Index active_rows)
{
    solver_type solver;
    REQUIRE(solver.setup(problem).has_value());

    const auto result = solver.solve(ctrlpp::nlp_update<double>{Eigen::VectorXd::Zero(2)});

    REQUIRE(result.status == ctrlpp::solve_status::optimal);
    REQUIRE(result.stop_criterion == ctrlpp::nlp_stop_criterion::stationarity);
    REQUIRE(result.x.size() == 2);
    CHECK((result.x - Eigen::Vector2d::Constant(0.5 * optimum_sum)).lpNorm<Eigen::Infinity>() <= optimum_tolerance(active_rows));
}

void check_solves_to(double lower, double upper, double optimum_sum, Eigen::Index active_rows)
{
    SECTION("analytic Jacobian")
    {
        check_path(make_problem(lower, upper, jacobian_path::analytic), optimum_sum, active_rows);
    }

    SECTION("finite-difference Jacobian")
    {
        check_path(make_problem(lower, upper, jacobian_path::finite_difference), optimum_sum, active_rows);
    }
}

}

TEST_CASE("an active upper side holds the row at its upper bound", "[argmin][inequality]")
{
    check_solves_to(-inf, 1.0, 1.0, 1);
}

TEST_CASE("an active lower side holds the row at its lower bound", "[argmin][inequality]")
{
    check_solves_to(3.0, inf, 3.0, 1);
}

TEST_CASE("a two-sided row holds at whichever side is active", "[argmin][inequality]")
{
    SECTION("the upper side")
    {
        check_solves_to(-1.0, 1.0, 1.0, 1);
    }

    SECTION("the lower side")
    {
        check_solves_to(3.0, 5.0, 3.0, 1);
    }
}

TEST_CASE("an inactive row leaves the unconstrained optimum", "[argmin][inequality]")
{
    SECTION("upper side only")
    {
        check_solves_to(-inf, 5.0, 2.0, 0);
    }

    SECTION("lower side only")
    {
        check_solves_to(-5.0, inf, 2.0, 0);
    }

    SECTION("both sides")
    {
        check_solves_to(-5.0, 5.0, 2.0, 0);
    }
}
