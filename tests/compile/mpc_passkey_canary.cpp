#include "ctrlpp/mpc.h"

enum class canary_setup_error
{
    rejected,
};

struct canary_qp_solver
{
    using scalar_type = double;

    auto setup(const ctrlpp::qp_problem<double>&) -> ctrlpp::expected<void, canary_setup_error> { return {}; }

    auto solve(const ctrlpp::qp_update<double>&) -> ctrlpp::qp_result<double> { return {}; }
};

using canary_controller = ctrlpp::mpc<double, 1, 1, canary_qp_solver>;

int main()
{
    ctrlpp::discrete_state_space<double, 1, 1, 1> system{};
    system.A.setIdentity();
    system.B.setOnes();
    system.C.setIdentity();
    system.D.setZero();
    const ctrlpp::mpc_config<double, 1, 1> cfg{};
#if defined(CTRLPP_CANARY_MPC_DIRECT_CONSTRUCTION)
    canary_controller controller{{}, system, cfg, ctrlpp::Matrix<double, 1, 1>::Identity(), ctrlpp::detail::mpc_terminal_weights<double, 1, 1>{}, canary_qp_solver{}};
    return controller.trajectory().has_value() ? 1 : 0;
#else
    return canary_controller::create(system, cfg).has_value() ? 0 : 1;
#endif
}
