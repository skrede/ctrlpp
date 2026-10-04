#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_ARM_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_ARM_H

#include "instance.h"

#include <Eigen/Dense>

#include <string>
#include <cstdint>

namespace ctrlpp::bench::qp_solvers
{

enum class solver_claim : std::uint8_t
{
    solved,
    iteration_cap,
    infeasible,
    error
};

/// x in the instance's variable order; y in OSQP's sign convention (positive on an active upper
/// bound) over the instance's row order, empty when the solver exposes no multipliers.
struct arm_answer
{
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    std::int32_t iterations;
    std::string status;
    solver_claim claim;
    bool fallback;
};

/// setup() poses the step-0 tick and may solve it; update() hands over the next tick's vectors through
/// the solver's own update call, or stages them into the arm's buffers when the solve call takes them;
/// solve() solves from the previous solution. answer() is read outside every timed span.
class qp_arm
{
public:
    virtual ~qp_arm() = default;

    virtual void setup(const qp_instance& instance, const qp_tick& tick) = 0;

    virtual void update(const qp_tick& tick) = 0;

    virtual void solve() = 0;

    virtual void answer(arm_answer& out) = 0;
};

}

#endif
