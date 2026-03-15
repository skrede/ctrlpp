#include "ctrlpp/linalg_policy.h"
#include "ctrlpp/solver_policy.h"
#include "ctrlpp/observer_policy.h"

#include "naive_linalg.h"

// LinalgPolicy conformance
static_assert(ctrlpp::LinalgPolicy<NaiveLinalg>);
static_assert(!ctrlpp::LinalgPolicy<int>);

// SolverPolicy conformance
struct NullSolver {
    using solver_tag = void;
};

static_assert(ctrlpp::SolverPolicy<NullSolver>);
static_assert(!ctrlpp::SolverPolicy<int>);

// ObserverPolicy conformance
static_assert(ctrlpp::ObserverPolicy<ctrlpp::NullObserver>);
static_assert(!ctrlpp::ObserverPolicy<int>);

int main() { return 0; }
