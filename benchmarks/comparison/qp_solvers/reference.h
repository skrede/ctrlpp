#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_REFERENCE_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_REFERENCE_H

#include "cases.h"
#include "registry.h"
#include "hpipm_arm.h"

#include "qp/qp_reference.h"

#include <vector>
#include <cstdint>

namespace ctrlpp::bench::qp_solvers
{

/// The reference answer's inequality-row classes on one tick; found is false when no reference
/// solver reached the reference accuracy there.
struct reference_step
{
    bool found;
    std::vector<problems::qp::row_state> rows;
};

/// The forms one problem is posed in for its reference solves: the dense form (condensed, or the
/// dense-only problem) and, for the MPC cells, the sparse form HPIPM takes.
struct reference_forms
{
    const problem_case& dense;
    const problem_case* sparse;
};

/// Solves every tick of one sequence with the reference solvers in their fixed order (qpOASES in its
/// reliable preset on the dense form, HPIPM in its balanced mode on the sparse form, the active-set
/// solver on the dense form), keeping the first answer within the reference accuracy.
auto reference_sequence(const reference_forms& forms, regime kind, std::int32_t seed, std::int32_t steps,
                        hpipm_dual_signs signs) -> std::vector<reference_step>;

auto row_classes(const inequality_rows& rows) -> std::vector<problems::qp::row_state>;

/// -1 when the step has no reference answer.
auto binding_fraction(const reference_step& step) -> double;

/// -1 unless both steps have a reference answer.
auto churn_fraction(const reference_step& previous, const reference_step& current) -> double;

}

#endif
