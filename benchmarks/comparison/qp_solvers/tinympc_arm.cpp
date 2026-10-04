#include "tinympc_arm.h"
#include "tinympc_session.h"

#include <vector>
#include <optional>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

constexpr double tinympc_rho = 1.0;

auto tinympc_options_of(tier level) -> tinympc_options
{
    return {solver_tolerance(level), iteration_cap, tinympc_rho, false};
}

auto as_values(const Eigen::MatrixXd& matrix) -> std::vector<double>
{
    return {matrix.data(), matrix.data() + matrix.size()};
}

/// TinyMPC takes the stage costs as diagonals; the instance's are diagonal.
auto problem_of(const qp_instance& instance) -> tinympc_problem
{
    const ocp_structure& ocp = *instance.ocp;
    const auto& cell = ocp.cell;
    const std::int32_t states = (cell.horizon + 1) * cell.nx;
    const Eigen::VectorXd state_weight = instance.P.diagonal().head(cell.nx);
    const Eigen::VectorXd input_weight = instance.P.diagonal().segment(states, cell.nu);
    return {cell.nx,          cell.nu,         cell.horizon,          as_values(ocp.transition), as_values(ocp.actuation),
            as_values(state_weight), as_values(input_weight), as_values(ocp.target), ocp.input_limit, ocp.state_limit};
}

class tinympc_arm final : public qp_arm
{
public:
    explicit tinympc_arm(const tinympc_options& options)
        : x0_{}
        , states_(0)
        , inputs_(0)
        , status_(0)
        , options_(options)
        , session_{}
    {
    }

    void setup(const qp_instance& instance, const qp_tick& tick) override
    {
        session_.emplace(problem_of(instance), options_);
        const auto& cell = instance.ocp->cell;
        states_ = (cell.horizon + 1) * cell.nx;
        inputs_ = cell.horizon * cell.nu;
        x0_ = tick.x0;
    }

    void update(const qp_tick& tick) override
    {
        x0_ = tick.x0;
        session_->set_initial_state(x0_.data());
    }

    void solve() override
    {
        status_ = session_->solve();
    }

    /// TinyMPC's knot-major x and u are the instance's stacked states and inputs.
    void answer(arm_answer& out) override
    {
        out.x.resize(states_ + inputs_);
        session_->read(out.x.data(), out.x.data() + states_);
        out.y.resize(0);
        out.iterations = session_->iterations();
        out.status = status_ == 0 ? "TINY_SOLVED" : status_ == 1 ? "TINY_MAX_ITER" : "TINY_ERROR";
        out.claim = status_ == 0 ? solver_claim::solved : status_ == 1 ? solver_claim::iteration_cap : solver_claim::error;
        out.fallback = false;
    }

private:
    Eigen::VectorXd x0_;
    std::int32_t states_;
    std::int32_t inputs_;
    std::int32_t status_;
    tinympc_options options_;
    std::optional<tinympc_session> session_;
};

}

auto make_tinympc_arm(tier level) -> std::unique_ptr<qp_arm>
{
    return std::make_unique<tinympc_arm>(tinympc_options_of(level));
}

void describe_tinympc_arm(json_writer& out, tier level)
{
    describe_tinympc(out, tinympc_options_of(level));
}

auto tinympc_unit_eigen_version() -> std::string
{
    return tinympc_eigen_version();
}

}
