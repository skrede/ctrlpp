#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_DRIVER_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_DRIVER_H

#include "registry.h"

namespace ctrlpp::bench::qp_solvers
{

/// One process of one cohort: --describe writes the settings and registry the manifest is built from,
/// --binding-gate writes the gate over the full registry, and otherwise the convention probes run and,
/// if every one passes, one repetition's blocks are measured into <out>/<cohort>_rep<r>.csv.
auto run_driver(cohort group, int argc, char** argv) -> int;

}

#endif
