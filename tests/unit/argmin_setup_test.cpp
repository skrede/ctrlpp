#include "ctrlpp/nmpc.h"

#include "ctrlpp/mpc/nlp_types.h"
#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/argmin_solver.h"
#include "ctrlpp/mpc/argmin_policies.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>

#include <span>
#include <cmath>
#include <limits>
#include <cstddef>
#include <variant>

namespace
{

constexpr double inf = std::numeric_limits<double>::infinity();

using capped_solver  = ctrlpp::argmin_solver<double, ctrlpp::argmin_nw_sqp, true, 3, 3>;
using static_problem = ctrlpp::nlp_problem_static<double, 3>;

struct linear_rows
{
    Eigen::MatrixXd coefficients;
    Eigen::VectorXd lower;
    Eigen::VectorXd upper;
};

auto as_vector(std::span<const double> x) -> Eigen::Map<const Eigen::VectorXd>
{
    return {x.data(), static_cast<Eigen::Index>(x.size())};
}

// Minimizes 0.5 |x - target|^2 subject to lower <= A x <= upper, row by row.
template <typename Problem>
auto make_projection(const Eigen::VectorXd& target, const linear_rows& rows) -> Problem
{
    Problem problem{};
    problem.n_vars        = static_cast<int>(target.size());
    problem.n_constraints = static_cast<int>(rows.coefficients.rows());
    problem.cost_cb          = [target](std::span<const double> x) { return 0.5 * (as_vector(x) - target).squaredNorm(); };
    problem.gradient_cb      = [target](std::span<const double> x, std::span<double> g)
    { Eigen::Map<Eigen::VectorXd>(g.data(), target.size()) = as_vector(x) - target; };
    problem.constraints_cb = [a = rows.coefficients](std::span<const double> x, std::span<double> c)
    { Eigen::Map<Eigen::VectorXd>(c.data(), a.rows()) = a * as_vector(x); };
    // The raw Jacobian is column-major, n_constraints x n_vars.
    problem.constraint_jacobian_cb = [a = rows.coefficients](std::span<const double>, std::span<double> jac)
    { Eigen::Map<Eigen::MatrixXd>(jac.data(), a.rows(), a.cols()) = a; };
    problem.x_lower = decltype(problem.x_lower)::Constant(target.size(), -10.0);
    problem.x_upper = decltype(problem.x_upper)::Constant(target.size(), 10.0);
    problem.c_lower = rows.lower;
    problem.c_upper = rows.upper;
    return problem;
}

auto projection_target() -> Eigen::VectorXd
{
    return Eigen::Vector3d{1.0, 0.0, 0.0};
}

// One equality row and one two-sided row: three partitioned rows. The two-sided row's upper side is active at the
// optimum (0.75, 0.25, 0) with multiplier 0.25, and its lower side is not.
auto three_partitioned_rows() -> linear_rows
{
    Eigen::MatrixXd a(2, 3);
    a << 1.0, 1.0, 1.0, 1.0, -1.0, 0.0;
    return {a, Eigen::Vector2d{1.0, -0.5}, Eigen::Vector2d{1.0, 0.5}};
}

auto three_partitioned_rows_active() -> Eigen::MatrixXd
{
    return three_partitioned_rows().coefficients;
}

auto three_partitioned_rows_optimum() -> Eigen::Vector3d
{
    return {0.75, 0.25, 0.0};
}

// The three rows above plus the one-sided row x2 <= 1, inactive at the same optimum.
auto four_partitioned_rows() -> linear_rows
{
    Eigen::MatrixXd a(3, 3);
    a << 1.0, 1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 1.0;
    return {a, Eigen::Vector3d{1.0, -0.5, -inf}, Eigen::Vector3d{1.0, 0.5, 1.0}};
}

auto two_two_sided_rows() -> linear_rows
{
    Eigen::MatrixXd a(2, 3);
    a << 1.0, -1.0, 0.0, 0.0, 1.0, -1.0;
    return {a, Eigen::Vector2d{-0.5, -0.5}, Eigen::Vector2d{0.5, 0.5}};
}

// On the active set, x - x* = K^-1 r for the KKT matrix K of the identity Hessian and the active rows, and kkt_tol
// bounds every entry of the residual r, so |x - x*| <= sqrt(dim r) * kkt_tol / sigma_min(K).
auto optimum_tolerance(const Eigen::MatrixXd& active, double kkt_tol) -> double
{
    const Eigen::Index n = active.cols();
    const Eigen::Index m = active.rows();
    Eigen::MatrixXd kkt  = Eigen::MatrixXd::Zero(n + m, n + m);
    kkt.topLeftCorner(n, n).setIdentity();
    kkt.topRightCorner(n, m)   = active.transpose();
    kkt.bottomLeftCorner(m, n) = active;
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(kkt);
    return std::sqrt(static_cast<double>(n + m)) * kkt_tol / svd.singularValues().minCoeff();
}

template <typename Solver>
auto solve_from_origin(Solver& solver, Eigen::Index n) -> ctrlpp::nlp_result<double>
{
    return solver.solve(ctrlpp::nlp_update<double>{Eigen::VectorXd::Zero(n)});
}

void check_reaches(const ctrlpp::nlp_result<double>& result, const Eigen::VectorXd& optimum, const Eigen::MatrixXd& active)
{
    const double tolerance = optimum_tolerance(active, ctrlpp::argmin_settings<double>{}.kkt_tol);
    REQUIRE(result.status == ctrlpp::solve_status::optimal);
    REQUIRE(result.stop_criterion == ctrlpp::nlp_stop_criterion::stationarity);
    REQUIRE(result.x.size() == optimum.size());
    CHECK((result.x - optimum).lpNorm<Eigen::Infinity>() <= tolerance);
}


auto equality_rows(const Eigen::MatrixXd& a, const Eigen::VectorXd& value) -> linear_rows
{
    return {a, value, value};
}

constexpr std::size_t nx = 2;
constexpr std::size_t nu = 1;
constexpr std::size_t nh = 5;
constexpr int nv         = static_cast<int>((nh + 1) * nx + nh * nu);
constexpr int posed_rows = static_cast<int>(nx * (nh + 1));

struct double_integrator
{
    auto operator()(const ctrlpp::Vector<double, nx>& x, const ctrlpp::Vector<double, nu>& u) const -> ctrlpp::Vector<double, nx>
    {
        constexpr double dt = 0.1;
        return {x[0] + dt * x[1], x[1] + dt * u[0]};
    }
};

template <int MaxM>
using integrator_controller = ctrlpp::nmpc_static<double, nx, nu, nh, ctrlpp::argmin_solver<double, ctrlpp::argmin_nw_sqp, true, nv, MaxM>, double_integrator>;

auto bounded_input_config() -> ctrlpp::nmpc_config<double, nx, nu>
{
    ctrlpp::nmpc_config<double, nx, nu> config{};
    config.horizon = static_cast<int>(nh);
    config.u_min   = ctrlpp::Vector<double, nu>::Constant(-1.0);
    config.u_max   = ctrlpp::Vector<double, nu>::Constant(1.0);
    return config;
}
}

TEST_CASE("exactly MaxM partitioned rows set up and solve", "[argmin][setup][cap]")
{
    capped_solver solver;
    const auto problem = make_projection<static_problem>(projection_target(), three_partitioned_rows());

    REQUIRE(solver.setup(problem).has_value());
    check_reaches(solve_from_origin(solver, 3), three_partitioned_rows_optimum(), three_partitioned_rows_active());
}

TEST_CASE("one partitioned row past MaxM is refused", "[argmin][setup][cap]")
{
    capped_solver solver;
    const auto problem = make_projection<static_problem>(projection_target(), four_partitioned_rows());

    const auto refused = solver.setup(problem);

    REQUIRE_FALSE(refused.has_value());
    CHECK(refused.error() == ctrlpp::argmin_setup_error::constraint_cap_exceeded);
}

TEST_CASE("a two-sided row counts twice against MaxM", "[argmin][setup][cap]")
{
    capped_solver solver;
    const auto problem = make_projection<static_problem>(projection_target(), two_two_sided_rows());
    REQUIRE(problem.n_constraints == 2);

    const auto refused = solver.setup(problem);

    REQUIRE_FALSE(refused.has_value());
    CHECK(refused.error() == ctrlpp::argmin_setup_error::constraint_cap_exceeded);
}

TEST_CASE("a cap refusal leaves the solver as it was", "[argmin][setup][cap]")
{
    capped_solver solver;
    const auto fits     = make_projection<static_problem>(projection_target(), three_partitioned_rows());
    const auto one_over = make_projection<static_problem>(projection_target(), four_partitioned_rows());
    const auto two_over = make_projection<static_problem>(projection_target(), two_two_sided_rows());

    REQUIRE_FALSE(solver.setup(one_over).has_value());
    REQUIRE(solver.setup(fits).has_value());
    check_reaches(solve_from_origin(solver, 3), three_partitioned_rows_optimum(), three_partitioned_rows_active());

    SECTION("a refused setup keeps the problem that fits bound")
    {
        REQUIRE_FALSE(solver.setup(two_over).has_value());
        check_reaches(solve_from_origin(solver, 3), three_partitioned_rows_optimum(), three_partitioned_rows_active());
    }
}

TEST_CASE("a problem whose n_vars differs from the fixed NV is refused", "[argmin][setup][shape]")
{
    capped_solver solver;
    auto problem   = make_projection<static_problem>(projection_target(), three_partitioned_rows());
    problem.n_vars = 2;

    const auto refused = solver.setup(problem);

    REQUIRE_FALSE(refused.has_value());
    CHECK(refused.error() == ctrlpp::argmin_setup_error::dimension_mismatch);
}

TEST_CASE("more equality rows than decision variables are refused", "[argmin][setup][shape]")
{
    using dynamic_solver = ctrlpp::argmin_solver<double, ctrlpp::argmin_nw_sqp>;
    const Eigen::VectorXd target = Eigen::Vector2d{1.0, 0.0};

    SECTION("three equality rows on two variables")
    {
        Eigen::MatrixXd a(3, 2);
        a << 1.0, 0.0, 0.0, 1.0, 1.0, 1.0;
        const auto problem = make_projection<ctrlpp::nlp_problem<double>>(target, equality_rows(a, Eigen::Vector3d{0.5, 0.5, 1.0}));
        dynamic_solver solver;

        const auto refused = solver.setup(problem);

        REQUIRE_FALSE(refused.has_value());
        CHECK(refused.error() == ctrlpp::argmin_setup_error::equality_rows_exceed_dimension);
    }

    SECTION("two equality rows on two variables set up and solve")
    {
        const Eigen::MatrixXd a = Eigen::Matrix2d::Identity();
        const auto problem      = make_projection<ctrlpp::nlp_problem<double>>(target, equality_rows(a, Eigen::Vector2d{0.5, 0.5}));
        dynamic_solver solver;

        REQUIRE(solver.setup(problem).has_value());
        check_reaches(solve_from_origin(solver, 2), Eigen::Vector2d{0.5, 0.5}, a);
    }
}

TEST_CASE("setup refuses with the first failing check", "[argmin][setup][shape]")
{
    SECTION("constraint rows on the unconstrained bridge come before a dimension mismatch")
    {
        ctrlpp::argmin_solver<double, ctrlpp::argmin_lbfgsb, false, 3> solver;
        auto problem   = make_projection<static_problem>(projection_target(), three_partitioned_rows());
        problem.n_vars = 2;

        const auto refused = solver.setup(problem);

        REQUIRE_FALSE(refused.has_value());
        CHECK(refused.error() == ctrlpp::argmin_setup_error::constraints_unsupported);
    }

    SECTION("a dimension mismatch comes before the cap")
    {
        capped_solver solver;
        auto problem   = make_projection<static_problem>(projection_target(), four_partitioned_rows());
        problem.n_vars = 2;

        const auto refused = solver.setup(problem);

        REQUIRE_FALSE(refused.has_value());
        CHECK(refused.error() == ctrlpp::argmin_setup_error::dimension_mismatch);
    }

    SECTION("the cap comes before the equality-row count")
    {
        Eigen::MatrixXd a(4, 3);
        a << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0;
        capped_solver solver;
        const auto problem = make_projection<static_problem>(projection_target(), equality_rows(a, Eigen::Vector4d{0.25, 0.25, 0.5, 1.0}));

        const auto refused = solver.setup(problem);

        REQUIRE_FALSE(refused.has_value());
        CHECK(refused.error() == ctrlpp::argmin_setup_error::constraint_cap_exceeded);
    }
}

TEST_CASE("nmpc_static refuses a solver whose cap is below the posed rows", "[argmin][setup][controller]")
{
    static_assert(posed_rows == 12);

    SECTION("one row short of the posed count is refused at construction")
    {
        const auto created = integrator_controller<posed_rows - 1>::create(double_integrator{}, bounded_input_config());

        REQUIRE_FALSE(created.has_value());
        const auto* cause = std::get_if<ctrlpp::argmin_setup_error>(&created.error());
        REQUIRE(cause != nullptr);
        CHECK(*cause == ctrlpp::argmin_setup_error::constraint_cap_exceeded);
    }

    SECTION("the posed count is accepted and solves")
    {
        auto created = integrator_controller<posed_rows>::create(double_integrator{}, bounded_input_config());

        REQUIRE(created.has_value());
        CHECK(created->solve(ctrlpp::Vector<double, nx>{1.0, 0.0}).has_value());
    }
}
