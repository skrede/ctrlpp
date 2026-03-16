#include "ctrlpp/implot/app.h"
#include "ctrlpp/implot/sim_harness.h"
#include "ctrlpp/implot/signal_recorder.h"
#include "ctrlpp/implot/scrolling_plot.h"
#include "ctrlpp/implot/static_plot.h"

#include "ctrlpp/pid.h"
#include "ctrlpp/pid_policies.h"
#include "ctrlpp/eigen_linalg.h"

#include "imgui.h"

int main()
{
    using Policy = ctrlpp::EigenLinalgPolicy;
    using Vec = Eigen::Matrix<double, 1, 1>;

    // Outer loop: temperature PI with anti-windup and IAE
    using OuterPid = ctrlpp::Pid<double, 1, 1, 1, Policy,
        ctrlpp::AntiWindup<ctrlpp::BackCalc>,
        ctrlpp::PerfAssessment<ctrlpp::IAE>>;

    OuterPid::config_type outer_cfg{};
    outer_cfg.kp = Vec::Constant(1.5);
    outer_cfg.ki = Vec::Constant(0.3);
    outer_cfg.kd = Vec::Constant(0.0);
    outer_cfg.output_min = Vec::Constant(0.0);
    outer_cfg.output_max = Vec::Constant(100.0);

    OuterPid outer(outer_cfg);

    // Inner loop: heater power PID with derivative filter
    using InnerPid = ctrlpp::Pid<double, 1, 1, 1, Policy,
        ctrlpp::DerivFilter>;

    InnerPid::config_type inner_cfg{};
    inner_cfg.kp = Vec::Constant(3.0);
    inner_cfg.ki = Vec::Constant(1.0);
    inner_cfg.kd = Vec::Constant(0.2);
    inner_cfg.output_min = Vec::Constant(0.0);
    inner_cfg.output_max = Vec::Constant(100.0);
    inner_cfg.template policy<ctrlpp::DerivFilter>().n = {10.0};

    InnerPid inner(inner_cfg);

    // Plant states
    double temperature = 20.0;
    double heater_power = 0.0;
    double temp_setpoint = 60.0;
    constexpr double dt = 0.01;

    ctrlpp::implot::SignalRecorder recorder(6000);

    ctrlpp::implot::SimHarness sim(dt,
        [&](double t) {
            Vec sp_outer = Vec::Constant(temp_setpoint);
            Vec meas_outer = Vec::Constant(temperature);
            Vec inner_output = Vec::Constant(heater_power);

            auto heater_sp_vec = outer.compute(sp_outer, meas_outer, dt, inner_output);
            double heater_sp = heater_sp_vec[0];

            Vec sp_inner = Vec::Constant(heater_sp);
            Vec meas_inner = Vec::Constant(heater_power);
            auto control = inner.compute(sp_inner, meas_inner, dt);

            heater_power = 0.8 * heater_power + 0.2 * control[0];
            temperature = 0.99 * temperature + 0.01 * heater_power;

            recorder.record("temp_setpoint", t, temp_setpoint);
            recorder.record("temperature", t, temperature);
            recorder.record("heater_sp", t, heater_sp);
            recorder.record("heater_power", t, heater_power);
            recorder.record("control", t, control[0]);
        },
        [&] {
            outer.reset();
            inner.reset();
            temperature = 20.0;
            heater_power = 0.0;
            recorder.clear();
        });

    ctrlpp::implot::App app(1400, 800, "Advanced PID -- Cascade Temperature Control");

    app.run([&] {
        sim.advance(ImGui::GetIO().DeltaTime);

        ImGui::Begin("Simulation");
        sim.draw_controls();
        ImGui::End();

        // Outer loop parameters
        ImGui::Begin("Outer Loop (Temperature PI)");

        auto sp_f = static_cast<float>(temp_setpoint);
        if (ImGui::SliderFloat("Setpoint", &sp_f, 0.0f, 100.0f))
            temp_setpoint = sp_f;

        ImGui::Separator();

        auto okp_f = static_cast<float>(outer_cfg.kp[0]);
        auto oki_f = static_cast<float>(outer_cfg.ki[0]);
        bool outer_changed = false;
        outer_changed |= ImGui::SliderFloat("Kp##outer", &okp_f, 0.0f, 10.0f);
        outer_changed |= ImGui::SliderFloat("Ki##outer", &oki_f, 0.0f, 5.0f);
        if (outer_changed) {
            outer_cfg.kp[0] = okp_f;
            outer_cfg.ki[0] = oki_f;
            outer.set_params(outer_cfg);
        }

        ImGui::Separator();
        auto iae = outer.metric<ctrlpp::IAE>();
        ImGui::Text("IAE: %.4f", iae[0]);

        ImGui::End();

        // Inner loop parameters
        ImGui::Begin("Inner Loop (Heater PID)");

        auto ikp_f = static_cast<float>(inner_cfg.kp[0]);
        auto iki_f = static_cast<float>(inner_cfg.ki[0]);
        auto ikd_f = static_cast<float>(inner_cfg.kd[0]);
        bool inner_changed = false;
        inner_changed |= ImGui::SliderFloat("Kp##inner", &ikp_f, 0.0f, 20.0f);
        inner_changed |= ImGui::SliderFloat("Ki##inner", &iki_f, 0.0f, 10.0f);
        inner_changed |= ImGui::SliderFloat("Kd##inner", &ikd_f, 0.0f, 5.0f);

        auto n_val = static_cast<float>(
            inner_cfg.template policy<ctrlpp::DerivFilter>().n[0]);
        inner_changed |= ImGui::SliderFloat("DerivFilter N", &n_val, 1.0f, 50.0f);

        auto imax_f = static_cast<float>(inner_cfg.output_max[0]);
        inner_changed |= ImGui::SliderFloat("Output Max##inner", &imax_f, 0.0f, 200.0f);

        if (inner_changed) {
            inner_cfg.kp[0] = ikp_f;
            inner_cfg.ki[0] = iki_f;
            inner_cfg.kd[0] = ikd_f;
            inner_cfg.template policy<ctrlpp::DerivFilter>().n[0] = n_val;
            inner_cfg.output_max[0] = imax_f;
            inner.set_params(inner_cfg);
        }

        ImGui::End();

        // SVG export
        ImGui::Begin("Export");
        if (ImGui::Button("Export SVG"))
            ctrlpp::implot::write_svg("advanced_pid.svg", recorder,
                {"temp_setpoint", "temperature", "heater_sp", "heater_power", "control"});
        ImGui::End();

        // Plots
        auto t = static_cast<float>(sim.sim_time());
        ctrlpp::implot::scrolling_plot(
            {"Temperature", "Time (s)", "deg", 20.0f},
            recorder, {"temp_setpoint", "temperature"}, t);
        ctrlpp::implot::scrolling_plot(
            {"Heater Power", "Time (s)", "W", 20.0f},
            recorder, {"heater_sp", "heater_power"}, t);
        ctrlpp::implot::scrolling_plot(
            {"Control Signal", "Time (s)", "u", 20.0f},
            recorder, {"control"}, t);
    });
}
