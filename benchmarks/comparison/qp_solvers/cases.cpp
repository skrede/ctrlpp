#include "cases.h"

#include "qp/mpc_condensing.h"
#include "qp/terminal_riccati.h"
#include "qp/mpc_initial_states.h"

namespace ctrlpp::bench::qp_solvers
{

namespace
{

namespace qp = problems::qp;

auto spec_of(const qp::mpc_sweep_cell& cell, const qp::mpc_variant& variant) -> qp::mpc_tracking_spec
{
    return {cell, variant, qp::terminal_weight(cell, qp::sweep_plant(cell))};
}

auto ocp_of(const qp::mpc_tracking_spec& spec) -> ocp_structure
{
    const qp::mpc_sweep_plant plant = qp::sweep_plant(spec.cell);
    return {spec.cell, plant.transition, plant.actuation, spec.terminal_weight, qp::tracking_target(spec),
            spec.variant.input_limit, spec.variant.state_limit};
}

auto dense_instance(const Eigen::MatrixXd& P, const Eigen::MatrixXd& A, std::int32_t box_rows) -> qp_instance
{
    return {P, A, P.sparseView(), A.sparseView(), box_rows, std::nullopt};
}

class sparse_case final : public problem_case
{
public:
    explicit sparse_case(const qp::mpc_tracking_spec& spec)
        : spec_(spec)
        , problem_(qp::build_tracking_problem(spec))
        , instance_{Eigen::MatrixXd(problem_.P), Eigen::MatrixXd(problem_.A), problem_.P, problem_.A, 0, ocp_of(spec)}
    {
    }

    auto instance() const -> const qp_instance& override
    {
        return instance_;
    }

    void fill(const Eigen::VectorXd& x0, qp_tick& tick) const override
    {
        tick.q = problem_.q;
        tick.l = problem_.l;
        tick.u = problem_.u;
        qp::set_initial_state(x0, tick.l, tick.u);
        tick.x0 = x0;
    }

    auto rows_of(const qp_tick& tick, const Eigen::VectorXd& x) const -> inequality_rows override
    {
        const Eigen::Index first = problem_.layout.n_dyn;
        const Eigen::Index count = problem_.A.rows() - first;
        return {(problem_.A * x).tail(count), tick.l.tail(count), tick.u.tail(count)};
    }

    auto state_dimension() const -> std::int32_t override
    {
        return spec_.cell.nx;
    }

    auto amplitude() const -> double override
    {
        return spec_.variant.x0_amplitude;
    }

private:
    qp::mpc_tracking_spec spec_;
    qp::mpc_tracking_problem problem_;
    qp_instance instance_;
};

class condensed_case final : public problem_case
{
public:
    explicit condensed_case(const qp::mpc_tracking_spec& spec)
        : spec_(spec)
        , posed_(qp::condense_tracking(spec))
        , instance_(dense_instance(posed_.P, posed_.A, spec.cell.horizon * spec.cell.nu))
    {
    }

    auto instance() const -> const qp_instance& override
    {
        return instance_;
    }

    void fill(const Eigen::VectorXd& x0, qp_tick& tick) const override
    {
        tick.q.resize(posed_.P.rows());
        tick.l = posed_.l;
        tick.u = posed_.u;
        qp::condense_initial_state(spec_, posed_, x0, tick.q, tick.l, tick.u);
        tick.x0 = x0;
    }

    /// The state rows are shifted back by the free response Phi x0 that condensing moved into their bounds.
    auto rows_of(const qp_tick& tick, const Eigen::VectorXd& x) const -> inequality_rows override
    {
        Eigen::VectorXd shift = Eigen::VectorXd::Zero(tick.l.size());
        const std::int32_t state_rows = qp::tracking_state_rows(spec_);
        if(state_rows > 0)
            shift.tail(state_rows) = posed_.actuated_phi * tick.x0;
        return {posed_.A * x + shift, tick.l + shift, tick.u + shift};
    }

    auto state_dimension() const -> std::int32_t override
    {
        return spec_.cell.nx;
    }

    auto amplitude() const -> double override
    {
        return spec_.variant.x0_amplitude;
    }

private:
    qp::mpc_tracking_spec spec_;
    qp::condensed_tracking posed_;
    qp_instance instance_;
};

}

auto make_sparse_case(const qp::mpc_sweep_cell& cell, const qp::mpc_variant& variant) -> std::unique_ptr<problem_case>
{
    return std::make_unique<sparse_case>(spec_of(cell, variant));
}

auto make_condensed_case(const qp::mpc_sweep_cell& cell, const qp::mpc_variant& variant) -> std::unique_ptr<problem_case>
{
    return std::make_unique<condensed_case>(spec_of(cell, variant));
}

auto initial_state(const problem_case& source, regime kind, std::int32_t seed, std::int32_t step) -> Eigen::VectorXd
{
    const std::int32_t nx = source.state_dimension();
    switch(kind)
    {
        case regime::fixed:
            return qp::fixed_initial_state(nx, source.amplitude());
        case regime::closed_loop:
            return qp::closed_loop_initial_state(nx, source.amplitude(), step);
        default:
            return qp::churning_initial_state(nx, source.amplitude(), seed, step);
    }
}

}
