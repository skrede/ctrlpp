#include "telemetry.h"
#include "family_report.h"

#include <cstdio>

namespace ctrlpp {

void report_family(const family_result &result)
{
    char body[kLineBytes];
    std::snprintf(body, sizeof(body), "family %s value=%.17g golden=%.17g departure=%.3e bound=%.3e verdict=%s allocations=unmeasured", result.name, result.value, result.golden,
                  result.departure, result.bound, result.pass ? "PASS" : "FAIL");
    send_line(stream_tag::verdict, body);
}

void report_refusal(const char *name, const char *reason)
{
    char body[kLineBytes];
    std::snprintf(body, sizeof(body), "family %s verdict=REFUSED -- %s", name, reason);
    send_line(stream_tag::verdict, body);
}

void report_unmeasured_allocations()
{
    send_line(stream_tag::verdict, "alloc steady_state allocations=unmeasured -- this image carries no allocation sensor");
}

void report_absent_predictive()
{
    send_line(stream_tag::verdict, "family predictive verdict=ABSENT -- the image was built without the predictive backend");
}

void report_golden(const control_loop_demo<float> &demo, const golden_verdict &verdict)
{
    const double final_norm = static_cast<double>(demo.x.norm());
    char body[kLineBytes];
    std::snprintf(body, sizeof(body), "gain check %s: dev = [%.3e, %.3e], tol = %.3e", verdict.gain_pass ? "PASS" : "FAIL", verdict.gain_departure[0], verdict.gain_departure[1],
                  verdict.gain_tolerance);
    send_line(stream_tag::verdict, body);
    std::snprintf(body, sizeof(body), "settling: final |x| = %.6e (host = %.6e, err = %.3e, tol = %.3e)", final_norm, kHostFinalNorm, verdict.norm_departure, verdict.norm_tolerance);
    send_line(stream_tag::verdict, body);
    report_family({"control", final_norm, kHostFinalNorm, verdict.norm_departure, verdict.norm_tolerance, verdict.pass});
    std::snprintf(body, sizeof(body), "golden diff %s (gain %s, final norm %s)", verdict.pass ? "PASS" : "FAIL", verdict.gain_pass ? "PASS" : "FAIL",
                  verdict.norm_pass ? "PASS" : "FAIL");
    send_line(stream_tag::verdict, body);
}

}
