#include "ctrlpp/nmpc.h"
#include "ctrlpp/types.h"

#include "ctrlpp/mpc/nlp_types.h"
#include "ctrlpp/mpc/diagnostics.h"
#include "ctrlpp/mpc/argmin_solver.h"
#include "ctrlpp/mpc/argmin_policies.h"

#include <argmin/solver/convergence.h>

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

#include <cstddef>
#include <cstdint>

namespace
{

constexpr std::size_t nx = 2;
constexpr std::size_t nu = 1;
constexpr std::size_t nh = 5;
constexpr int nv         = static_cast<int>((nh + 1) * nx + nh * nu);
constexpr int max_m      = static_cast<int>(nx * (nh + 1));
constexpr double dt      = 0.1;
constexpr int run_length = 221;

struct double_integrator
{
    ctrlpp::Vector<double, nx> operator()(const ctrlpp::Vector<double, nx>& x, const ctrlpp::Vector<double, nu>& u) const
    {
        return ctrlpp::Vector<double, nx>{x[0] + dt * x[1], x[1] + dt * u[0]};
    }
};

using solver_type     = ctrlpp::argmin_solver<double, ctrlpp::argmin_nw_sqp, true, nv, max_m>;
using controller_type = ctrlpp::nmpc_static<double, nx, nu, nh, solver_type, double_integrator>;

ctrlpp::nmpc_config<double, nx, nu> make_config()
{
    ctrlpp::nmpc_config<double, nx, nu> config;
    config.horizon = static_cast<int>(nh);
    config.Q       = 10.0 * Eigen::Matrix2d::Identity();
    config.R       = 0.1 * Eigen::Matrix<double, 1, 1>::Identity();
    return config;
}

struct run_record
{
    Eigen::Vector2d final_state;
    std::int64_t iterations;
    int stationarity_stops;
    int iteration_limit_stops;
    bool all_solved;
};

run_record run(const ctrlpp::argmin_settings<double>& settings)
{
    controller_type controller{double_integrator{}, make_config(), solver_type{settings}};
    run_record record{Eigen::Vector2d{1.0, 0.0}, 0, 0, 0, true};
    for(int k = 0; k < run_length; ++k)
    {
        const auto solved = controller.solve(record.final_state);
        if(!solved.has_value())
        {
            record.all_solved = false;
            return record;
        }
        const auto diagnostics = controller.diagnostics();
        record.iterations += diagnostics.iterations;
        record.stationarity_stops += diagnostics.stop_criterion == ctrlpp::nlp_stop_criterion::stationarity ? 1 : 0;
        record.iteration_limit_stops += diagnostics.stop_criterion == ctrlpp::nlp_stop_criterion::iteration_limit ? 1 : 0;
        record.final_state = double_integrator{}(record.final_state, solved->input);
    }
    return record;
}

}

TEST_CASE("kkt_tol defaults to the backend's own stationarity threshold", "[argmin][settings]")
{
    REQUIRE(ctrlpp::argmin_settings<double>{}.kkt_tol == argmin::gradient_tolerance_criterion{}.threshold);
}

TEST_CASE("leaving kkt_tol at its default reproduces the run that sets the backend's threshold bit for bit", "[argmin][settings]")
{
    ctrlpp::argmin_settings<double> explicit_default;
    explicit_default.kkt_tol = argmin::gradient_tolerance_criterion{}.threshold;

    const run_record defaulted = run(ctrlpp::argmin_settings<double>{});
    const run_record pinned    = run(explicit_default);

    REQUIRE(defaulted.all_solved);
    REQUIRE(pinned.all_solved);
    REQUIRE(defaulted.iterations == pinned.iterations);
    REQUIRE(defaulted.final_state[0] == pinned.final_state[0]);
    REQUIRE(defaulted.final_state[1] == pinned.final_state[1]);
}

TEST_CASE("a tighter kkt_tol reaches the backend and buys more iterations", "[argmin][settings]")
{
    ctrlpp::argmin_settings<double> tight;
    tight.kkt_tol = 1e-9;

    const run_record defaulted = run(ctrlpp::argmin_settings<double>{});
    const run_record tightened = run(tight);

    REQUIRE(defaulted.all_solved);
    REQUIRE(tightened.all_solved);
    REQUIRE(tightened.iterations > defaulted.iterations);
}

TEST_CASE("the stopping test is passed through to the controller's diagnostics", "[argmin][diagnostics]")
{
    SECTION("every solve of the default run stops on stationarity")
    {
        const run_record defaulted = run(ctrlpp::argmin_settings<double>{});
        REQUIRE(defaulted.all_solved);
        REQUIRE(defaulted.stationarity_stops == run_length);
    }

    SECTION("a one-evaluation budget under an unreachable threshold stops on the budget")
    {
        ctrlpp::argmin_settings<double> starved;
        starved.kkt_tol  = 0.0;
        starved.max_eval = 1;
        const run_record limited = run(starved);
        REQUIRE(limited.all_solved);
        REQUIRE(limited.iteration_limit_stops > 0);
        REQUIRE(limited.stationarity_stops == 0);
    }

    SECTION("a QP-backed controller's diagnostics leave it unreported")
    {
        REQUIRE(ctrlpp::mpc_diagnostics<double>{}.stop_criterion == ctrlpp::nlp_stop_criterion::unreported);
    }
}
