#include "hpipm_arm.h"
#include "hpipm_stages.h"
#include "hpipm_settings.h"

#include <vector>
#include <cstddef>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

/// HPIPM's create functions place their structures inside caller-provided memory of a size they report.
struct hpipm_memory
{
    std::vector<std::byte> dim;
    std::vector<std::byte> qp;
    std::vector<std::byte> sol;
    std::vector<std::byte> arg;
    std::vector<std::byte> ws;
};

class hpipm_arm final : public qp_arm
{
public:
    hpipm_arm(const hpipm_config& config, hpipm_dual_signs signs)
        : status_(SUCCESS)
        , qp_{}
        , iterations_(0)
        , dim_{}
        , sol_{}
        , ws_{}
        , x0_{}
        , config_(config)
        , stages_{}
        , memory_{}
        , arg_{}
        , signs_(signs)
    {
    }

    void setup(const qp_instance& instance, const qp_tick& tick) override;

    void update(const qp_tick& tick) override
    {
        x0_ = tick.x0;
        d_ocp_qp_set_lbx(0, x0_.data(), &qp_);
        d_ocp_qp_set_ubx(0, x0_.data(), &qp_);
    }

    void solve() override
    {
        d_ocp_qp_ipm_solve(&qp_, &sol_, &arg_, &ws_);
        d_ocp_qp_ipm_get_status(&ws_, &status_);
        d_ocp_qp_ipm_get_iter(&ws_, &iterations_);
    }

    void answer(arm_answer& out) override;

private:
    int status_;
    d_ocp_qp qp_;
    int iterations_;
    d_ocp_qp_dim dim_;
    d_ocp_qp_sol sol_;
    d_ocp_qp_ipm_ws ws_;
    Eigen::VectorXd x0_;
    hpipm_config config_;
    hpipm_stages stages_;
    hpipm_memory memory_;
    d_ocp_qp_ipm_arg arg_;
    hpipm_dual_signs signs_;

    void create_dimensions();

    void read_primal(Eigen::VectorXd& x);

    void read_equality_duals(Eigen::VectorXd& y);

    void read_bound_duals(Eigen::VectorXd& y);
};

void hpipm_arm::create_dimensions()
{
    memory_.dim.resize(d_ocp_qp_dim_memsize(stages_.horizon));
    d_ocp_qp_dim_create(stages_.horizon, &dim_, memory_.dim.data());
    d_ocp_qp_dim_set_all(stages_.nx.data(), stages_.nu.data(), stages_.nbx.data(), stages_.nbu.data(), stages_.ng.data(),
                         stages_.ns.data(), &dim_);
    d_ocp_qp_dim_set_nbxe(0, stages_.nx.front(), &dim_);
}

void hpipm_arm::setup(const qp_instance& instance, const qp_tick& tick)
{
    stages_ = stages_of(*instance.ocp);
    create_dimensions();
    memory_.qp.resize(d_ocp_qp_memsize(&dim_));
    d_ocp_qp_create(&dim_, &qp_, memory_.qp.data());
    pose_stages(instance, tick, stages_, qp_);
    memory_.sol.resize(d_ocp_qp_sol_memsize(&dim_));
    d_ocp_qp_sol_create(&dim_, &sol_, memory_.sol.data());
    memory_.arg.resize(d_ocp_qp_ipm_arg_memsize(&dim_));
    d_ocp_qp_ipm_arg_create(&dim_, &arg_, memory_.arg.data());
    configure(config_, arg_);
    memory_.ws.resize(d_ocp_qp_ipm_ws_memsize(&dim_, &arg_));
    d_ocp_qp_ipm_ws_create(&dim_, &arg_, &ws_, memory_.ws.data());
    x0_ = tick.x0;
}

void hpipm_arm::read_primal(Eigen::VectorXd& x)
{
    const int N = stages_.horizon;
    const int nx = stages_.nx.front();
    const int nu = stages_.nu.front();
    x.resize((N + 1) * nx + N * nu);
    for(int k = 0; k <= N; ++k)
        d_ocp_qp_sol_get_x(k, &sol_, x.data() + k * nx);
    for(int k = 0; k < N; ++k)
        d_ocp_qp_sol_get_u(k, &sol_, x.data() + (N + 1) * nx + k * nu);
}

/// The stage-0 equality bound and the dynamics multipliers take the probe-settled signs.
void hpipm_arm::read_equality_duals(Eigen::VectorXd& y)
{
    const int N = stages_.horizon;
    const int nx = stages_.nx.front();
    Eigen::VectorXd lower(nx);
    Eigen::VectorXd upper(nx);
    d_ocp_qp_sol_get_lam_lbx(0, &sol_, lower.data());
    d_ocp_qp_sol_get_lam_ubx(0, &sol_, upper.data());
    y.head(nx) = signs_.initial_state * (upper - lower);
    for(int k = 0; k < N; ++k)
        d_ocp_qp_sol_get_pi(k, &sol_, y.data() + nx + k * nx);
    y.segment(nx, N * nx) *= signs_.dynamics;
}

/// Bound multipliers map as upper minus lower.
void hpipm_arm::read_bound_duals(Eigen::VectorXd& y)
{
    const int N = stages_.horizon;
    const int nu = stages_.nu.front();
    const int states = stages_.nbx.back();
    const int offset = (N + 1) * stages_.nx.front();
    Eigen::VectorXd lower(stages_.nx.front());
    Eigen::VectorXd upper(stages_.nx.front());
    for(int k = 0; k < N; ++k)
    {
        d_ocp_qp_sol_get_lam_lbu(k, &sol_, lower.data());
        d_ocp_qp_sol_get_lam_ubu(k, &sol_, upper.data());
        y.segment(offset + k * nu, nu) = upper.head(nu) - lower.head(nu);
    }
    for(int k = 1; k <= N && states > 0; ++k)
    {
        d_ocp_qp_sol_get_lam_lbx(k, &sol_, lower.data());
        d_ocp_qp_sol_get_lam_ubx(k, &sol_, upper.data());
        y.segment(offset + N * nu + (k - 1) * states, states) = upper.head(states) - lower.head(states);
    }
}

void hpipm_arm::answer(arm_answer& out)
{
    read_primal(out.x);
    out.y.resize(out.x.size() + stages_.horizon * stages_.nbx.back());
    read_equality_duals(out.y);
    read_bound_duals(out.y);
    out.iterations = iterations_;
    out.status = hpipm_status_token(status_);
    out.claim = hpipm_claim(status_);
    out.fallback = false;
}

}

auto make_hpipm_arm(tier level, hpipm_dual_signs signs) -> std::unique_ptr<qp_arm>
{
    return std::make_unique<hpipm_arm>(config_of(level), signs);
}

auto make_hpipm_reference(hpipm_dual_signs signs) -> std::unique_ptr<qp_arm>
{
    return std::make_unique<hpipm_arm>(reference_config(), signs);
}

}
