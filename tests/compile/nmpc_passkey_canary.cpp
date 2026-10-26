#include "ctrlpp/nmpc.h"

#include <Eigen/Core>

enum class canary_setup_error
{
    rejected,
};

struct canary_solver
{
    using scalar_type = double;

    auto setup(const ctrlpp::nlp_problem<double>&) -> ctrlpp::expected<void, canary_setup_error> { return {}; }

    auto solve(const ctrlpp::nlp_update<double>&) -> ctrlpp::nlp_result<double>
    {
        return ctrlpp::nlp_result<double>{.status = ctrlpp::solve_status::optimal, .x = Eigen::VectorXd::Zero(3), .objective = 0.0, .solve_time = 0.0, .iterations = 0, .primal_residual = 0.0};
    }
};

struct canary_dynamics
{
    auto operator()(const ctrlpp::Vector<double, 1>& x, const ctrlpp::Vector<double, 1>&) const -> ctrlpp::Vector<double, 1> { return x; }
};

using canary_controller = ctrlpp::nmpc_dynamic<double, 1, 1, canary_solver, canary_dynamics>;

int main()
{
    const ctrlpp::nmpc_config<double, 1, 1> cfg{};
#ifdef CTRLPP_CANARY_NMPC_DYNAMIC_DIRECT_CONSTRUCTION
    canary_controller controller{{}, canary_dynamics{}, cfg, canary_solver{}};
    return controller.problem().n_vars > 0 ? 0 : 1;
#else
    const auto controller = canary_controller::create(canary_dynamics{}, cfg);
    return controller.has_value() ? 0 : 1;
#endif
}
