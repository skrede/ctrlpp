#include "hpipm_stages.h"

#include <numeric>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

struct sparse_rows
{
    int inputs;
    int states;
};

void pose_costs(const qp_instance& instance, const qp_tick& tick, int k, int N, d_ocp_qp& qp)
{
    const auto& cell = instance.ocp->cell;
    const int nx = cell.nx;
    const int nu = cell.nu;
    const int states = (N + 1) * nx;
    Eigen::MatrixXd Q = instance.P.block(k * nx, k * nx, nx, nx);
    Eigen::VectorXd q = tick.q.segment(k * nx, nx);
    d_ocp_qp_set_Q(k, Q.data(), &qp);
    d_ocp_qp_set_q(k, q.data(), &qp);
    if(k == N)
        return;
    Eigen::MatrixXd R = instance.P.block(states + k * nu, states + k * nu, nu, nu);
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(nu, nx);
    Eigen::VectorXd r = tick.q.segment(states + k * nu, nu);
    d_ocp_qp_set_R(k, R.data(), &qp);
    d_ocp_qp_set_S(k, S.data(), &qp);
    d_ocp_qp_set_r(k, r.data(), &qp);
}

void pose_dynamics(const ocp_structure& ocp, int k, d_ocp_qp& qp)
{
    Eigen::MatrixXd A = ocp.transition;
    Eigen::MatrixXd B = ocp.actuation;
    Eigen::VectorXd b = Eigen::VectorXd::Zero(ocp.cell.nx);
    d_ocp_qp_set_A(k, A.data(), &qp);
    d_ocp_qp_set_B(k, B.data(), &qp);
    d_ocp_qp_set_b(k, b.data(), &qp);
}

void pose_bounds(const qp_tick& tick, const sparse_rows& rows, int k, hpipm_stages& stages, d_ocp_qp& qp)
{
    const int nu = stages.nu[0];
    if(k < stages.horizon)
    {
        Eigen::VectorXd lower = tick.l.segment(rows.inputs + k * nu, nu);
        Eigen::VectorXd upper = tick.u.segment(rows.inputs + k * nu, nu);
        d_ocp_qp_set_idxbu(k, stages.leading_inputs.data(), &qp);
        d_ocp_qp_set_lbu(k, lower.data(), &qp);
        d_ocp_qp_set_ubu(k, upper.data(), &qp);
    }
    if(k == 0 || stages.nbx[static_cast<std::size_t>(k)] == 0)
        return;
    Eigen::VectorXd lower = tick.l.segment(rows.states + (k - 1) * nu, nu);
    Eigen::VectorXd upper = tick.u.segment(rows.states + (k - 1) * nu, nu);
    d_ocp_qp_set_idxbx(k, stages.leading_inputs.data(), &qp);
    d_ocp_qp_set_lbx(k, lower.data(), &qp);
    d_ocp_qp_set_ubx(k, upper.data(), &qp);
}

void pose_initial_state(const qp_tick& tick, hpipm_stages& stages, d_ocp_qp& qp)
{
    Eigen::VectorXd x0 = tick.x0;
    d_ocp_qp_set_idxbx(0, stages.leading_states.data(), &qp);
    d_ocp_qp_set_idxbxe(0, stages.leading_states.data(), &qp);
    d_ocp_qp_set_lbx(0, x0.data(), &qp);
    d_ocp_qp_set_ubx(0, x0.data(), &qp);
}

}

auto stages_of(const ocp_structure& ocp) -> hpipm_stages
{
    const int N = ocp.cell.horizon;
    const auto stages = static_cast<std::size_t>(N + 1);
    hpipm_stages out{N, std::vector<int>(stages, ocp.cell.nx), std::vector<int>(stages, ocp.cell.nu),
                     std::vector<int>(stages, ocp.state_limit ? ocp.cell.nu : 0), std::vector<int>(stages, ocp.cell.nu),
                     std::vector<int>(stages, 0), std::vector<int>(stages, 0), std::vector<int>(static_cast<std::size_t>(ocp.cell.nx)),
                     std::vector<int>(static_cast<std::size_t>(ocp.cell.nu))};
    out.nu.back() = 0;
    out.nbu.back() = 0;
    out.nbx.front() = ocp.cell.nx;
    std::iota(out.leading_states.begin(), out.leading_states.end(), 0);
    std::iota(out.leading_inputs.begin(), out.leading_inputs.end(), 0);
    return out;
}

void pose_stages(const qp_instance& instance, const qp_tick& tick, hpipm_stages& stages, d_ocp_qp& qp)
{
    const ocp_structure& ocp = *instance.ocp;
    const int N = stages.horizon;
    const int inputs = (N + 1) * ocp.cell.nx;
    const sparse_rows rows{inputs, inputs + N * ocp.cell.nu};
    for(int k = 0; k <= N; ++k)
    {
        if(k < N)
            pose_dynamics(ocp, k, qp);
        pose_costs(instance, tick, k, N, qp);
        pose_bounds(tick, rows, k, stages, qp);
    }
    pose_initial_state(tick, stages, qp);
}

}
