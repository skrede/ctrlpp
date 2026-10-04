#include "probe_problems.h"

#include "qp/mpc_tracking_problem.h"

#include <ranges>
#include <stdexcept>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

namespace qp = problems::qp;

/// The probe is only a probe if its exact answer has an upper-active row, a lower-active row and an
/// equality row.
void require_pattern(const probe_problem& probe)
{
    const auto& rows = probe.exact.rows;
    const bool upper = std::ranges::count(rows, qp::row_state::at_upper) > 0;
    const bool lower = std::ranges::count(rows, qp::row_state::at_lower) > 0;
    const bool equality = (probe.tick.l.array() == probe.tick.u.array()).any();
    if(probe.exact.x.size() == 0 || !upper || !lower || !equality)
        throw std::logic_error("convention probe does not pose an active upper bound, an active lower bound and an equality row");
}

auto finish(qp_instance instance, qp_tick tick) -> probe_problem
{
    probe_problem probe{std::move(instance), std::move(tick), {}};
    probe.exact = qp::solve_by_enumeration(program_of(probe.instance, probe.tick));
    require_pattern(probe);
    return probe;
}

}

auto dense_probe() -> probe_problem
{
    Eigen::MatrixXd P = Eigen::MatrixXd::Identity(3, 3);
    Eigen::MatrixXd A(4, 3);
    A << 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1;
    qp_instance instance{P, A, P.sparseView(), A.sparseView(), 3, std::nullopt};
    qp_tick tick{Eigen::Vector3d(-3.0, 3.0, 1.0), Eigen::Vector4d(-1.0, -1.0, -1.0, 0.0), Eigen::Vector4d(1.0, 1.0, 1.0, 0.0), {}};
    return finish(std::move(instance), std::move(tick));
}

auto sparse_probe() -> probe_problem
{
    const qp::mpc_sweep_cell cell{1, 1, 2};
    const qp::mpc_variant variant{"probe", 1.0, std::nullopt, 0.0, 0.0};
    const qp::mpc_tracking_spec spec{cell, variant, Eigen::MatrixXd::Identity(1, 1)};
    const qp::mpc_tracking_problem problem = qp::build_tracking_problem(spec);
    const qp::mpc_sweep_plant plant = qp::sweep_plant(cell);
    const ocp_structure ocp{cell, plant.transition, plant.actuation, spec.terminal_weight, Eigen::VectorXd::Zero(1), 1.0, std::nullopt};
    qp_instance instance{Eigen::MatrixXd(problem.P), Eigen::MatrixXd(problem.A), problem.P, problem.A, 0, ocp};
    qp_tick tick{problem.q, problem.l, problem.u, Eigen::VectorXd::Constant(1, 0.5)};
    tick.q(1) = -3.0;
    tick.q(2) = 0.5;
    qp::set_initial_state(tick.x0, tick.l, tick.u);
    return finish(std::move(instance), std::move(tick));
}

}
