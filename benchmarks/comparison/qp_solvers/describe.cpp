#include "describe.h"
#include "arm_factory.h"
#include "tinympc_arm.h"
#include "alloc_sensor.h"
#include "instruction_counter.h"

#include "qp/terminal_riccati.h"
#include "qp/mpc_tracking_cells.h"

#include <Eigen/Core>

#include <format>
#include <limits>
#include <fstream>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

namespace qp = problems::qp;

auto registered_once(const arm_spec& arm) -> bool
{
    return arm.group == cohort::dense || arm.kind != arm_kind::osqp || arm.forced_polish;
}

void write_settings(json_writer& out)
{
    out.open("arm_settings");
    for(const arm_spec& arm : arm_registry)
    {
        if(!registered_once(arm))
            continue;
        out.open(arm.id);
        for(const tier level : tiers)
            if(runs_at(arm, level))
            {
                out.open(name_of(level));
                describe_arm(out, arm, level);
                out.close();
            }
        out.close();
    }
    out.close();
    out.open("update_api");
    for(const arm_spec& arm : arm_registry)
        if(registered_once(arm))
            out.text(arm.id, name_of(arm.update));
    out.close();
}

void write_expressibility(json_writer& out)
{
    out.open("expressible");
    for(const arm_spec& arm : arm_registry)
    {
        out.open(std::format("{}/{}", name_of(arm.group), arm.id));
        for(const qp::mpc_sweep_cell& cell : qp::tracking_cells())
            out.flag(qp::cell_id(cell), arm.kind != arm_kind::tinympc || qp::is_long_horizon(cell));
        if(arm.group == cohort::dense)
            out.flag("dense_mpc_shaped", true);
        out.close();
    }
    out.close();
}

void write_variant(json_writer& out, const qp::mpc_variant& variant)
{
    out.open(std::format("{}_variant", variant.name));
    out.number("input_limit", variant.input_limit);
    out.number("state_limit", variant.state_limit.value_or(std::numeric_limits<double>::quiet_NaN()));
    out.number("tracking_target", variant.tracking_target);
    out.number("x0_amplitude", variant.x0_amplitude);
    out.close();
}

void write_tolerances(json_writer& out)
{
    out.open("tiers");
    out.number("mpc", mpc_threshold);
    out.number("high", high_threshold);
    out.close();
    out.number("convention_probe_tolerance", probe_tolerance);
    out.number("reference_solution_kkt_max", reference_kkt_max);
    out.number("binding_row_tolerance", binding_row_tolerance);
    out.number("binding_median_fraction_min", binding_median_fraction_min);
    out.number("binding_state_row_step_fraction_min", binding_state_row_step_fraction_min);
    out.integer("osqp_adaptive_rho_interval", osqp_adaptive_rho_interval);
}

void write_constants(json_writer& out)
{
    out.open("constants");
    out.integer("repetitions", repetitions);
    out.integer("steps_per_sequence", steps_per_sequence);
    out.integer("fixed_warmup_resolves", fixed_warmup_resolves);
    out.integer("fixed_timed_resolves", steps_per_sequence);
    write_tolerances(out);
    write_variant(out, qp::nominal_variant());
    write_variant(out, qp::tight_variant());
    out.number("legacy_terminal_weight", qp::legacy_terminal_weight);
    out.open("terminal_riccati");
    out.number("rho", qp::terminal_riccati.rho);
    out.integer("max_iterations", qp::terminal_riccati.max_iterations);
    out.number("gain_tolerance", qp::terminal_riccati.gain_tolerance);
    out.close();
    out.close();
}

/// Each entry in shortest round-trip form, so the reader recovers every binary64 exactly.
void write_terminal_weights(json_writer& out)
{
    out.open("terminal_weights");
    for(const qp::mpc_sweep_cell& cell : qp::tracking_cells())
    {
        if(!qp::is_long_horizon(cell))
            continue;
        const Eigen::MatrixXd weight = qp::terminal_weight(cell, qp::sweep_plant(cell));
        out.open(qp::cell_id(cell));
        for(Eigen::Index i = 0; i < weight.size(); ++i)
            out.number(std::format("{}", i), weight.data()[i]);
        out.number("asymmetry", (weight - weight.transpose()).cwiseAbs().maxCoeff());
        out.close();
    }
    out.close();
}

void write_environment(json_writer& out)
{
    out.open("environment");
    out.text("eigen_harness", std::format("{}.{}.{}", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION));
    out.text("eigen_tinympc_unit", tinympc_unit_eigen_version());
    out.flag("alloc_sensor", alloc_sensor_present());
    out.flag("instruction_counter_rdpmc", instruction_counter().available());
    out.close();
}

}

auto write_description(const std::string& path) -> bool
{
    std::ofstream file(path);
    json_writer out(file);
    out.open();
    write_settings(out);
    write_expressibility(out);
    write_constants(out);
    write_terminal_weights(out);
    write_environment(out);
    out.close();
    return file.good();
}

}
