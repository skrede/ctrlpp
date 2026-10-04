#include "golden_verdict.h"
#include "predictive_demo.h"
#include "golden_reference.h"
#include "predictive_bound.h"
#include "predictive_verdict.h"
#include "predictive_reference.h"

#include <cmath>
#include <cstdio>
#include <cstddef>
#include <cstdlib>

namespace {

using namespace ctrlpp;

void print_constants(const predictive_constants &c)
{
    std::printf("kHostPredictiveK0            = %.17g\n", c.k0);
    std::printf("kHostPredictiveK1            = %.17g\n", c.k1);
    std::printf("kHostPredictiveInputRow      = %.17g\n", c.input_row);
    std::printf("kHostPredictiveKktCondition  = %.17g\n", c.kkt_condition);
    std::printf("kHostPredictiveJacobianSigma = %.17g\n", c.jacobian_sigma);
    std::printf("kHostPredictiveJacobianNorm  = %.17g\n", c.jacobian_norm);
    std::printf("kHostPredictiveCost          = %.17g\n", c.cost);
}

// Two double runs of the same computation may part by twice what one may part
// from the exact value; a header constant further than that from this run is
// stale.
bool header_current(const predictive_constants &c, predictive_workspace &workspace)
{
    const double gain  = 2.0 * predictive_gain_departure();
    const double row   = 2.0 * (predictive_input_row() - kHostPredictiveInputRow);
    const double sigma = 2.0 * (kHostPredictiveJacobianSigma - predictive_jacobian_sigma());
    const double cost  = 2.0 * host_cost_departure(workspace);
    const bool current = within_bound(std::abs(c.k0 - kHostPredictiveK0), gain) && within_bound(std::abs(c.k1 - kHostPredictiveK1), gain) &&
            within_bound(std::abs(c.input_row - kHostPredictiveInputRow), row) && within_bound(std::abs(c.jacobian_sigma - kHostPredictiveJacobianSigma), sigma) &&
            within_bound(std::abs(c.cost - kHostPredictiveCost), cost);
    const double relative = (predictive_input_row() - kHostPredictiveInputRow) / kHostPredictiveInputRow;
    const double norm     = 2.0 * counted_departure_bound<double>(static_cast<double>(kPredictiveNv), c.jacobian_norm);
    const bool indicators = within_bound(std::abs(c.kkt_condition - kHostPredictiveKktCondition), 2.0 * relative * c.kkt_condition) &&
            within_bound(std::abs(c.jacobian_norm - kHostPredictiveJacobianNorm), norm);
    std::printf("predictive header constants %s\n", current && indicators ? "agree with this run" : "are STALE");
    return current && indicators;
}

bool report_host_run(predictive_workspace &workspace)
{
    predictive_demo<double> demo;
    predictive_record record{};
    const auto ran = run_predictive(demo, record, [](std::size_t) {}, [](std::size_t) {});
    if(!ran.has_value())
    {
        std::fprintf(stderr, "predictive solve refused on the host: %s\n", describe(ran.error()));
        return false;
    }
    const predictive_verdict v = judge_predictive(record, workspace);
    std::printf("predictive host run: stationarity stops %u of %zu, worst gate %.3e\n", record.stationarity_stops, kPredictiveRunSteps, v.worst_gate);
    std::printf("predictive host cost = %.17g, departure = %.3e, bound = %.3e (%.3e of the cost), host share = %.3e\n", v.cost, v.departure, v.bound, v.bound / kHostPredictiveCost,
                host_cost_departure(workspace));
    std::printf("predictive host verdict = %s\n", v.pass ? "PASS" : "FAIL");
    return v.pass;
}

}

int main()
{
    const predictive_constants constants = compute_predictive_constants();
    print_constants(constants);
    predictive_workspace workspace = make_predictive_workspace();
    const bool current             = header_current(constants, workspace);
    const bool host                = report_host_run(workspace);
    return current && host ? EXIT_SUCCESS : EXIT_FAILURE;
}
