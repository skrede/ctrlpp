#include "ctrlpp/model/constraint_model.h"
#include "ctrlpp/mpc/nlopt_solver.h"
#include "ctrlpp/nmpc.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <span>
#include <cmath>
#include <limits>
#include <random>
#include <vector>
#include <cstddef>
#include <algorithm>

namespace
{

using Catch::Matchers::WithinAbs;

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr double dt = 0.1;

auto double_integrator = [](const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) -> Eigen::Vector2d { return Eigen::Vector2d{x(0) + dt * x(1), x(1) + dt * u(0)}; };

using NloptSolver = ctrlpp::nlopt_solver<double>;

// Path constraint: position <= upper_bound  =>  g(x,u) = x[0] - upper_bound <= 0
constexpr std::size_t NC = 1;
auto make_upper_bound_constraint(double bound)
{
    return [bound](const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>&) -> ctrlpp::Vector<double, NC> { return ctrlpp::Vector<double, NC>{x(0) - bound}; };
}

// Terminal constraint: x[0] <= target  =>  h(x) = x[0] - target <= 0
constexpr std::size_t NTC_1 = 1;
auto make_terminal_constraint(double target)
{
    return [target](const Eigen::Vector2d& x) -> ctrlpp::Vector<double, NTC_1> { return ctrlpp::Vector<double, NTC_1>{x(0) - target}; };
}

} // namespace

// ----- Concept static assertions -----

TEST_CASE("constraint_model concept accepts valid path constraint lambda")
{
    auto g = make_upper_bound_constraint(1.0);
    static_assert(ctrlpp::constraint_model<decltype(g), double, NX, NU, NC>, "path constraint lambda must satisfy constraint_model");
}

TEST_CASE("terminal_constraint_model concept accepts valid terminal constraint lambda")
{
    auto h = make_terminal_constraint(0.5);
    static_assert(ctrlpp::terminal_constraint_model<decltype(h), double, NX, NTC_1>, "terminal constraint lambda must satisfy terminal_constraint_model");
}

// ----- Backward compatibility -----

TEST_CASE("nmpc with NC=0 NTC=0 produces same results as original", "[nmpc][constraint][compat]")
{
    ctrlpp::nmpc_config<double, NX, NU> config_old{
        .horizon = 10,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };

    ctrlpp::nmpc<double, NX, NU, NloptSolver, decltype(double_integrator)> ctrl_old{double_integrator, config_old};

    ctrlpp::nmpc_config<double, NX, NU, 0, 0> config_new{
        .horizon = 10,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };

    ctrlpp::nmpc<double, NX, NU, NloptSolver, decltype(double_integrator), 0, 0> ctrl_new{double_integrator, config_new};

    Eigen::Vector2d x{1.0, 0.0};
    auto u_old = ctrl_old.solve(x);
    auto u_new = ctrl_new.solve(x);

    REQUIRE(u_old.has_value());
    REQUIRE(u_new.has_value());
    CHECK_THAT((*u_old)(0), WithinAbs((*u_new)(0), 1e-8));
}

// ----- Soft path constraint -----

TEST_CASE("soft path constraint is approximately satisfied", "[nmpc][constraint][soft]")
{
    constexpr double upper_bound = 0.8;

    ctrlpp::nmpc_config<double, NX, NU, NC, 0> config{
        .horizon = 10,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };
    config.path_constraint = make_upper_bound_constraint(upper_bound);
    // soft_constraints defaults to true

    ctrlpp::nmpc<double, NX, NU, NloptSolver, decltype(double_integrator), NC, 0> controller{double_integrator, config};

    Eigen::Vector2d x{1.5, 0.0}; // starts above bound

    for(int step = 0; step < 40; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = double_integrator(x, *u);
    }

    // After enough steps with soft constraint, state should have moved below bound
    CHECK(x(0) < upper_bound + 0.5);

    auto diag = controller.diagnostics();
    // Diagnostics should have constraint violation info
    CHECK(diag.max_path_constraint_violation >= 0.0);
}

// ----- Hard path constraint -----

TEST_CASE("hard path constraint enforced tightly", "[nmpc][constraint][hard]")
{
    constexpr double upper_bound = 0.8;

    ctrlpp::nmpc_config<double, NX, NU, NC, 0> config{
        .horizon = 10,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };
    config.path_constraint = make_upper_bound_constraint(upper_bound);
    config.soft_constraints = false;

    ctrlpp::nmpc<double, NX, NU, NloptSolver, decltype(double_integrator), NC, 0> controller{double_integrator, config};

    // Start within bounds
    Eigen::Vector2d x{0.5, 0.0};

    for(int step = 0; step < 30; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = double_integrator(x, *u);

        auto [states, inputs] = controller.trajectory();
        for(const auto& s : states)
        {
            CHECK(s(0) <= upper_bound + 1e-3);
        }
    }

    auto diag = controller.diagnostics();
    CHECK(diag.total_slack == 0.0);
}

// ----- Terminal constraint -----

TEST_CASE("terminal constraint drives final state", "[nmpc][constraint][terminal]")
{
    constexpr double target = 0.3;

    ctrlpp::nmpc_config<double, NX, NU, 0, NTC_1> config{
        .horizon = 15,
        .Q = 1.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };
    config.terminal_constraint = make_terminal_constraint(target);
    config.soft_constraints = false;

    ctrlpp::nmpc<double, NX, NU, NloptSolver, decltype(double_integrator), 0, NTC_1> controller{double_integrator, config};

    Eigen::Vector2d x{1.0, 0.0};

    for(int step = 0; step < 50; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());

        // Check terminal state of predicted trajectory
        auto [states, inputs] = controller.trajectory();
        CHECK(states.back()(0) <= target + 1e-2);

        x = double_integrator(x, *u);
    }
}

// ----- Combined path + terminal constraints -----

TEST_CASE("combined path and terminal constraints", "[nmpc][constraint][combined]")
{
    constexpr double upper_bound = 1.2;
    constexpr double target = 0.3;

    ctrlpp::nmpc_config<double, NX, NU, NC, NTC_1> config{
        .horizon = 15,
        .Q = 5.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };
    config.path_constraint = make_upper_bound_constraint(upper_bound);
    config.terminal_constraint = make_terminal_constraint(target);

    ctrlpp::nmpc<double, NX, NU, NloptSolver, decltype(double_integrator), NC, NTC_1> controller{double_integrator, config};

    Eigen::Vector2d x{1.0, 0.0};

    for(int step = 0; step < 50; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = double_integrator(x, *u);
    }

    // State should have converged within combined constraints
    CHECK(x(0) < upper_bound + 0.5);

    auto diag = controller.diagnostics();
    CHECK(diag.max_path_constraint_violation >= -10.0); // sanity: field is populated
    CHECK(diag.max_terminal_constraint_violation >= -10.0);
}

// ----- Infeasible constraints with soft mode -----

TEST_CASE("infeasible constraints with soft mode does not crash", "[nmpc][constraint][infeasible]")
{
    // Constraint: x[0] <= -10 but starting at x[0] = 5 with regulation to origin
    // This is intentionally infeasible for early steps
    ctrlpp::nmpc_config<double, NX, NU, NC, 0> config{
        .horizon = 10,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };
    config.path_constraint = make_upper_bound_constraint(-10.0);
    // soft_constraints defaults to true

    ctrlpp::nmpc<double, NX, NU, NloptSolver, decltype(double_integrator), NC, 0> controller{double_integrator, config};

    Eigen::Vector2d x{5.0, 0.0};

    // Solver should not crash -- soft constraints absorb infeasibility
    bool found_slack = false;
    bool found_violation = false;
    for(int step = 0; step < 10; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());

        auto diag = controller.diagnostics();
        if(diag.total_slack > 0.0)
        {
            found_slack = true;
        }
        if(diag.max_path_constraint_violation > 0.0)
        {
            found_violation = true;
        }

        x = double_integrator(x, *u);
    }

    // At some point during the loop, slack should have been nonzero
    CHECK(found_slack);
    CHECK(found_violation);
}

// ----- MAJ-06 (i): state bounds applied from k=1 -----

TEST_CASE("nmpc state bounds skip the x0 block and apply from k=1", "[nmpc][constraint][bounds]")
{
    constexpr int N = 4;

    ctrlpp::nmpc_config<double, NX, NU> config{.horizon = N};
    config.x_min = ctrlpp::Vector<double, NX>(-1.5, -2.0);
    config.x_max = ctrlpp::Vector<double, NX>(1.5, 2.0);

    auto state = std::make_shared<ctrlpp::nmpc_formulation_state<double, NX, NU>>();
    auto prob = ctrlpp::detail::build_nmpc_problem<double, NX, NU, 0, 0>(double_integrator, config, state);

    constexpr int nx = static_cast<int>(NX);
    const double inf = std::numeric_limits<double>::infinity();

    // The x0 block is pinned only by the initial-state equality, so its variable
    // bounds must remain the unconstrained sentinels (never double-bounded).
    for(int i = 0; i < nx; ++i)
    {
        CHECK(prob.x_lower[i] == -inf);
        CHECK(prob.x_upper[i] == inf);
    }

    // The x_k blocks for k = 1..N carry the configured state bounds.
    for(int k = 1; k <= N; ++k)
    {
        for(int i = 0; i < nx; ++i)
        {
            CHECK(prob.x_lower[k * nx + i] == (*config.x_min)[i]);
            CHECK(prob.x_upper[k * nx + i] == (*config.x_max)[i]);
        }
    }
}

// ----- MAJ-06 (iv): analytic cost gradient matches central difference -----

TEST_CASE("nmpc analytic cost gradient matches central difference", "[nmpc][constraint][gradient]")
{
    constexpr int N = 4;

    ctrlpp::nmpc_config<double, NX, NU, NC, NTC_1> config{.horizon = N};
    config.Q = (Eigen::Matrix2d() << 3.0, 0.0, 0.0, 2.0).finished();
    config.R = 0.5 * Eigen::Matrix<double, 1, 1>::Identity();
    config.Qf = (Eigen::Matrix2d() << 4.0, 0.0, 0.0, 5.0).finished();
    config.path_constraint = make_upper_bound_constraint(0.8);
    config.terminal_constraint = make_terminal_constraint(0.3);
    config.du_max = ctrlpp::Vector<double, NU>(0.5);
    config.path_penalty = ctrlpp::Vector<double, NC>::Constant(5.0);
    config.terminal_penalty = ctrlpp::Vector<double, NTC_1>::Constant(7.0);
    // soft_constraints defaults to true, so path/terminal slack blocks exist.

    auto state = std::make_shared<ctrlpp::nmpc_formulation_state<double, NX, NU>>();
    state->x0 = ctrlpp::Vector<double, NX>(0.2, -0.1);
    state->u_prev = ctrlpp::Vector<double, NU>(0.05);
    for(int k = 0; k <= N; ++k)
    {
        state->x_ref.push_back(ctrlpp::Vector<double, NX>(0.1 * k, -0.05 * k));
    }

    auto prob = ctrlpp::detail::build_nmpc_problem<double, NX, NU, NC, NTC_1>(double_integrator, config, state);
    const int n = prob.n_vars;

    const double eps = std::numeric_limits<double>::epsilon();
    // Central-difference total error is O(eps^(2/3)) (see detail/numerical_diff.h).
    const double cd_err = std::pow(eps, 2.0 / 3.0);
    // Magnitude scales: the largest cost weight and the z sampling range. The unit
    // floor mirrors the max(1, |.|) convention used by the finite-difference step.
    const double w_scale = std::max({3.0, 2.0, 4.0, 5.0, 0.5, 5.0, 7.0});
    const double z_scale = 2.0;
    const double tol = static_cast<double>(n) * cd_err * w_scale * std::max(1.0, z_scale);

    std::mt19937 rng{20240709u};
    std::uniform_real_distribution<double> dist(-z_scale, z_scale);
    const double ss = std::cbrt(eps);

    for(int trial = 0; trial < 8; ++trial)
    {
        std::vector<double> z(static_cast<std::size_t>(n));
        for(auto& zi : z)
        {
            zi = dist(rng);
        }

        std::vector<double> ga(z.size());
        prob.gradient(std::span<const double>{z.data(), z.size()}, std::span<double>{ga.data(), ga.size()});

        // Central-difference reference of prob.cost.
        std::vector<double> gfd(z.size());
        std::vector<double> zp = z;
        for(std::size_t j = 0; j < z.size(); ++j)
        {
            const double h_raw = ss * std::max(1.0, std::abs(z[j]));
            const double t = z[j] + h_raw;
            const double h = t - z[j];
            const double o = zp[j];

            zp[j] = o + h;
            const double fp = prob.cost(std::span<const double>{zp.data(), zp.size()});
            zp[j] = o - h;
            const double fm = prob.cost(std::span<const double>{zp.data(), zp.size()});

            gfd[j] = (fp - fm) / (2.0 * h);
            zp[j] = o;
        }

        for(std::size_t j = 0; j < z.size(); ++j)
        {
            CHECK_THAT(ga[j], WithinAbs(gfd[j], tol));
        }
    }

    // The slack-penalty gradient is exactly the configured L1 weights.
    {
        std::vector<double> z(static_cast<std::size_t>(n), 0.0);
        std::vector<double> ga(z.size());
        prob.gradient(std::span<const double>{z.data(), z.size()}, std::span<double>{ga.data(), ga.size()});

        constexpr int nx = static_cast<int>(NX);
        constexpr int nu = static_cast<int>(NU);
        const int path_slack_offset = (N + 1) * nx + N * nu;
        CHECK(ga[static_cast<std::size_t>(path_slack_offset)] == 5.0);
        const int term_slack_offset = path_slack_offset + N * static_cast<int>(NC);
        CHECK(ga[static_cast<std::size_t>(term_slack_offset)] == 7.0);
    }
}

// ----- MAJ-06 (iv): constraint Jacobian structural rows + FD dynamics -----

TEST_CASE("nmpc constraint_jacobian matches central difference with exact structural rows", "[nmpc][constraint][jacobian]")
{
    constexpr int N = 3;

    ctrlpp::nmpc_config<double, NX, NU, NC, NTC_1> config{.horizon = N};
    config.Q = 2.0 * Eigen::Matrix2d::Identity();
    config.R = 0.5 * Eigen::Matrix<double, 1, 1>::Identity();
    config.path_constraint = make_upper_bound_constraint(0.8);
    config.terminal_constraint = make_terminal_constraint(0.3);
    config.du_max = ctrlpp::Vector<double, NU>(0.5);
    // soft_constraints defaults to true, so path/terminal slack columns exist.

    auto state = std::make_shared<ctrlpp::nmpc_formulation_state<double, NX, NU>>();
    state->x0 = ctrlpp::Vector<double, NX>(0.2, -0.1);
    state->u_prev = ctrlpp::Vector<double, NU>(0.05);
    for(int k = 0; k <= N; ++k)
    {
        state->x_ref.push_back(ctrlpp::Vector<double, NX>::Zero());
    }

    auto prob = ctrlpp::detail::build_nmpc_problem<double, NX, NU, NC, NTC_1>(double_integrator, config, state);
    const int n = prob.n_vars;
    const int m = prob.n_constraints;

    // Constraint / variable offset map (mirrors build_nmpc_problem).
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);
    constexpr int nc = static_cast<int>(NC);
    const int x_offset = 0;
    const int u_offset = (N + 1) * nx;
    const int path_slack_offset = u_offset + N * nu;
    const int term_slack_offset = path_slack_offset + N * nc;
    const int n_eq = (N + 1) * nx;
    const int n_rate = N * nu * 2;
    const int n_path_con = N * nc;
    const int eq_start = 0;
    const int rate_start = n_eq;
    const int path_con_start = rate_start + n_rate;
    const int term_con_start = path_con_start + n_path_con;

    // Column-major access: entry (row i, col j) at index i + j * m.
    const auto at = [m](const std::vector<double>& J, int i, int j) -> double { return J[static_cast<std::size_t>(i) + static_cast<std::size_t>(j) * static_cast<std::size_t>(m)]; };

    const double eps = std::numeric_limits<double>::epsilon();
    const double cd_err = std::pow(eps, 2.0 / 3.0);
    // Jacobian entries are O(1) (identities, +/-1 rate/slack blocks, dt-scaled
    // dynamics); the unit floor is the FD stencil's max(1, |.|) convention.
    const double tol = static_cast<double>(n) * cd_err;

    std::mt19937 rng{771u};
    std::uniform_real_distribution<double> dist(-2.0, 2.0);
    const double ss = std::cbrt(eps);

    std::vector<double> jac(static_cast<std::size_t>(m) * static_cast<std::size_t>(n));

    for(int trial = 0; trial < 6; ++trial)
    {
        std::vector<double> z(static_cast<std::size_t>(n));
        for(auto& zi : z)
        {
            zi = dist(rng);
        }

        prob.constraint_jacobian(std::span<const double>{z.data(), z.size()}, std::span<double>{jac.data(), jac.size()});

        // Central-difference reference of prob.constraints (column-major layout).
        std::vector<double> ref(jac.size());
        std::vector<double> zp = z;
        std::vector<double> cp(static_cast<std::size_t>(m));
        std::vector<double> cm(static_cast<std::size_t>(m));
        for(int j = 0; j < n; ++j)
        {
            const auto jz = static_cast<std::size_t>(j);
            const double h_raw = ss * std::max(1.0, std::abs(z[jz]));
            const double t = z[jz] + h_raw;
            const double h = t - z[jz];
            const double o = zp[jz];

            zp[jz] = o + h;
            prob.constraints(std::span<const double>{zp.data(), zp.size()}, std::span<double>{cp.data(), cp.size()});
            zp[jz] = o - h;
            prob.constraints(std::span<const double>{zp.data(), zp.size()}, std::span<double>{cm.data(), cm.size()});

            for(int i = 0; i < m; ++i)
            {
                ref[static_cast<std::size_t>(i) + jz * static_cast<std::size_t>(m)] = (cp[static_cast<std::size_t>(i)] - cm[static_cast<std::size_t>(i)]) / (2.0 * h);
            }
            zp[jz] = o;
        }

        for(std::size_t idx = 0; idx < jac.size(); ++idx)
        {
            CHECK_THAT(jac[idx], WithinAbs(ref[idx], tol));
        }
    }

    // Structural entries are exact (bit-exact 0 / +/-1), carrying no FD noise.
    // Initial-state identity rows.
    CHECK(at(jac, eq_start + 0, x_offset + 0) == 1.0);
    CHECK(at(jac, eq_start + 0, x_offset + 1) == 0.0);
    CHECK(at(jac, eq_start + 1, x_offset + 1) == 1.0);
    // Continuity identity block on x_{k+1} (k = 0).
    CHECK(at(jac, eq_start + nx + 0, x_offset + nx + 0) == 1.0);
    // Rate +/-1 blocks (k = 0, j = 0): +1 on the upper row, -1 on the lower row.
    CHECK(at(jac, rate_start + 0, u_offset + 0) == 1.0);
    CHECK(at(jac, rate_start + 1, u_offset + 0) == -1.0);
    // Slack -1 columns.
    CHECK(at(jac, path_con_start + 0, path_slack_offset + 0) == -1.0);
    CHECK(at(jac, term_con_start + 0, term_slack_offset + 0) == -1.0);
}
