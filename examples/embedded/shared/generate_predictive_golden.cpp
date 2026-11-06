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

namespace
{

using namespace ctrlpp;

void print_constants(const predictive_constants& c)
{
    std::printf("host_predictive_k0             = %.17g\n", c.k0);
    std::printf("host_predictive_k1             = %.17g\n", c.k1);
    std::printf("host_predictive_input_row      = %.17g\n", c.input_row);
    std::printf("host_predictive_kkt_condition  = %.17g\n", c.kkt_condition);
    std::printf("host_predictive_jacobian_sigma = %.17g\n", c.jacobian_sigma);
    std::printf("host_predictive_jacobian_norm  = %.17g\n", c.jacobian_norm);
    std::printf("host_predictive_cost           = %.17g\n", c.cost);
}

// Two double runs of the same computation may part by twice what one may part from the exact value; a
// header constant further than that from this run is stale.
bool header_current(const predictive_constants& c, predictive_workspace& workspace)
{
    const double gain = 2.0 * predictive_gain_departure();
    const double row = 2.0 * (predictive_input_row() - golden::host_predictive_input_row);
    const double sigma = 2.0 * (golden::host_predictive_jacobian_sigma - predictive_jacobian_sigma());
    const double cost = 2.0 * host_cost_departure(workspace);
    const bool current = within_bound(std::abs(c.k0 - golden::host_predictive_k0), gain) && within_bound(std::abs(c.k1 - golden::host_predictive_k1), gain) &&
            within_bound(std::abs(c.input_row - golden::host_predictive_input_row), row) && within_bound(std::abs(c.jacobian_sigma - golden::host_predictive_jacobian_sigma), sigma) &&
            within_bound(std::abs(c.cost - golden::host_predictive_cost), cost);
    const double relative = (predictive_input_row() - golden::host_predictive_input_row) / golden::host_predictive_input_row;
    const double norm = 2.0 * counted_departure_bound<double>(static_cast<double>(golden::predictive_nv), c.jacobian_norm);
    const bool indicators = within_bound(std::abs(c.kkt_condition - golden::host_predictive_kkt_condition), 2.0 * relative * c.kkt_condition) &&
            within_bound(std::abs(c.jacobian_norm - golden::host_predictive_jacobian_norm), norm);
    std::printf("predictive header constants %s\n", current && indicators ? "agree with this run" : "are STALE");
    return current && indicators;
}

bool report_host_run(predictive_workspace& workspace)
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
    std::printf("predictive host run: stationarity stops %u of %zu, worst gate %.3e\n", record.stationarity_stops, golden::predictive_run_steps, v.worst_gate);
    std::printf("predictive host cost = %.17g, departure = %.3e, bound = %.3e (%.3e of the cost), host share = %.3e\n", v.cost, v.departure, v.bound,
            v.bound / golden::host_predictive_cost, host_cost_departure(workspace));
    std::printf("predictive host verdict = %s\n", v.pass ? "PASS" : "FAIL");
    return v.pass;
}

}

int main()
{
    const predictive_constants constants = compute_predictive_constants();
    print_constants(constants);
    predictive_workspace workspace = make_predictive_workspace();
    const bool current = header_current(constants, workspace);
    const bool host = report_host_run(workspace);
    return current && host ? EXIT_SUCCESS : EXIT_FAILURE;
}
