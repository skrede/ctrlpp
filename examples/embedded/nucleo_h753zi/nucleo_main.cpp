// Armed before Eigen is parsed, or the Eigen allocation sentinel is never compiled in.
#include "eigen_alloc_sentinel.h"

// ctrlpp NUCLEO-H753ZI on-device leg (bare superloop, double FPU).
//
// CMSIS startup runs SystemInit + C++ static ctors, then calls main(). This
// prints the reset that started the run and the posture it is taken under, then
// either proves the allocation sensor is not blind (the evidence image) or proves
// the cycle counter live and times every family cold and warm (the timing image).
// It then designs the infinite-horizon LQR gain once, runs the closed loop with
// the sensor armed around each step, streams the trajectory's exact rows over
// USART3 -> ST-Link VCP, runs the estimation, signal-processing, trajectory and
// predictive families each inside its own armed window, diffs every family's
// double result against its host-double golden, and halts for the operator to
// read the console.

#include "posture.h"
#include "telemetry.h"
#include "timing_run.h"
#include "reset_cause.h"
#include "alloc_sensor.h"
#include "family_report.h"
#include "stack_watermark_board.h"

#include "dsp_demo.h"
#include "golden_verdict.h"
#include "estimation_demo.h"
#include "trajectory_demo.h"
#include "golden_reference.h"
#include "control_loop_demo.h"

#include <Eigen/Dense>

#include <cstdio>
#include <cstdint>

namespace ctrlpp {
void usart3_console_init() noexcept;
}

namespace {

using demo_type = ctrlpp::control_loop_demo<double>;

ctrlpp::golden_bound golden_workspace = ctrlpp::make_golden_bound();

void halt() noexcept
{
    for(;;)
    {
    }
}

// Runs before anything else executes a timed region, so the timing image's
// cold figures are first executions after reset.
void prove_instrument()
{
    if constexpr(ctrlpp::kAllocSensorPresent)
        ctrlpp::emit_allocation_canary(ctrlpp::run_alloc_canary());
    else
        ctrlpp::run_timing();
}

// The sensor is armed around step() alone, so the figure is the controller's
// and the telemetry formatting between steps is not charged to it.
void run_loop(demo_type &demo)
{
    ctrlpp::emit_header(ctrlpp::claim_class::exact);
    ctrlpp::alloc_sensor_reset();
    for(std::int32_t k = 0; k <= ctrlpp::kSteps; ++k)
    {
        ctrlpp::alloc_sensor_arm();
        const double u = demo.step();
        ctrlpp::alloc_sensor_disarm();
        ctrlpp::emit_exact_row(k, {k * ctrlpp::kDt, demo.x[0], demo.x[1], u});
    }
}

void report_golden(const demo_type &demo, const ctrlpp::golden_verdict &verdict, const ctrlpp::window_figures &figures)
{
    const double final_norm = demo.x.norm();
    ctrlpp::open_line(ctrlpp::stream_tag::verdict);
    std::printf("gain check %s: dev = [%.3e, %.3e], tol = %.3e\n", verdict.gain_pass ? "PASS" : "FAIL", verdict.gain_departure[0], verdict.gain_departure[1], verdict.gain_tolerance);
    ctrlpp::open_line(ctrlpp::stream_tag::verdict);
    std::printf("settling: final |x| = %.6e (host = %.6e, err = %.3e, tol = %.3e)\n", final_norm, ctrlpp::kHostFinalNorm, verdict.norm_departure, verdict.norm_tolerance);
    ctrlpp::report_family({"control", final_norm, ctrlpp::kHostFinalNorm, verdict.norm_departure, verdict.norm_tolerance, verdict.pass}, figures);
    ctrlpp::open_line(ctrlpp::stream_tag::verdict);
    std::printf("golden diff %s (gain %s, final norm %s)\n", verdict.pass ? "PASS" : "FAIL", verdict.gain_pass ? "PASS" : "FAIL", verdict.norm_pass ? "PASS" : "FAIL");
}

}

// The startup figure is read before main does anything else, so it covers reset
// and static construction; the program figure is read after the last family and
// the golden judgment, ahead of the closing report. The reset flags are cleared
// only once the posture line carrying them has been printed.
int main()
{
    const ctrlpp::stack_reading startup = ctrlpp::read_program_stack();
    ctrlpp::usart3_console_init();
    ctrlpp::emit_identity();
    ctrlpp::report_posture();
    ctrlpp::clear_reset_cause();
    ctrlpp::report_stack("startup", startup);
    ctrlpp::emit_header(ctrlpp::claim_class::spread);
    prove_instrument();

    auto demo = demo_type::make();
    if(!demo.has_value())
    {
        ctrlpp::open_line(ctrlpp::stream_tag::verdict);
        std::printf("family control verdict=REFUSED -- lqr_gain refused the plant on device: %s\n", ctrlpp::describe(demo.error()));
        halt();
    }

    ctrlpp::emit_setup_heap();
    ctrlpp::emit_gain(demo->K(0, 0), demo->K(0, 1));
    run_loop(*demo);
    const ctrlpp::window_figures control = ctrlpp::read_window();
    ctrlpp::report_allocations(control);
    ctrlpp::drive_family<ctrlpp::estimation_demo>("estimation", ctrlpp::kHostKalmanVelocityVariance, ctrlpp::estimation_tolerance<double>());
    ctrlpp::drive_family<ctrlpp::dsp_demo>("dsp", ctrlpp::kHostBiquadOutput, ctrlpp::dsp_tolerance<double>());
    ctrlpp::drive_family<ctrlpp::trajectory_demo>("trajectory", ctrlpp::kHostTrajectoryPosition, ctrlpp::trajectory_tolerance<double>());
    ctrlpp::drive_predictive_family();
    const ctrlpp::golden_verdict verdict = ctrlpp::judge_golden(*demo, golden_workspace);
    ctrlpp::report_stack("program", ctrlpp::read_program_stack());
    report_golden(*demo, verdict, control);
    halt();
}
