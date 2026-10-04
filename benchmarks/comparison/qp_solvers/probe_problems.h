#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_PROBE_PROBLEMS_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_PROBE_PROBLEMS_H

#include "instance.h"

#include "qp/qp_reference.h"

namespace ctrlpp::bench::qp_solvers
{

struct probe_problem
{
    qp_instance instance;
    qp_tick tick;
    problems::qp::reference_answer exact;
};

/// Three box rows on three variables and one equality row; the optimum holds the first variable on its
/// upper bound and the second on its lower bound, and x = 0 is feasible, as the active-set solver's
/// start requires. Every exact multiplier is of order one, the scale the absolute probe tolerance is
/// stated at.
auto dense_probe() -> probe_problem;

/// A one-state, one-input, two-step tracking problem in the sparse layout, whose linear state costs
/// put the first input on its upper bound and the second on its lower bound, with a nonzero initial
/// state so the initial-state and dynamics multipliers are exercised. Its multipliers are of order one.
auto sparse_probe() -> probe_problem;

}

#endif
