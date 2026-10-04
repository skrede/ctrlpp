#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_HPIPM_STAGES_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_HPIPM_STAGES_H

#include "instance.h"

extern "C"
{
#include <hpipm_d_ocp_qp.h>
#include <hpipm_d_ocp_qp_dim.h>
}

#include <vector>

namespace ctrlpp::bench::qp_solvers
{

/// Per-stage dimensions of the instance's OCP: the initial state as an equality bound on stage 0, the
/// input boxes on stages 0..N-1 and, with a state box, the actuated states on stages 1..N.
struct hpipm_stages
{
    int horizon;
    std::vector<int> nx;
    std::vector<int> nu;
    std::vector<int> nbx;
    std::vector<int> nbu;
    std::vector<int> ng;
    std::vector<int> ns;
    std::vector<int> leading_states;
    std::vector<int> leading_inputs;
};

auto stages_of(const ocp_structure& ocp) -> hpipm_stages;

void pose_stages(const qp_instance& instance, const qp_tick& tick, hpipm_stages& stages, d_ocp_qp& qp);

}

#endif
