#include "ctrlpp/implot/app.h"
#include "ctrlpp/implot/sim_harness.h"
#include "ctrlpp/implot/signal_recorder.h"
#include "ctrlpp/implot/scrolling_plot.h"
#include "ctrlpp/implot/static_plot.h"

#include "ctrlpp/pid.h"
#include "ctrlpp/eigen_linalg.h"

#include "imgui.h"

int main()
{
    using Policy = ctrlpp::EigenLinalgPolicy;
    using Pid = ctrlpp::Pid<double, 1, 1, 1, Policy>;
    using Vec = Pid::vector_t;

    // Controller config (mutable for slider binding)
    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);
    cfg.ki = Vec::Constant(1.0);
    cfg.kd = Vec::Constant(0.1);
    cfg.output_min = Vec::Constant(-10.0);
    cfg.output_max = Vec::Constant(10.0);

    Pid pid(cfg);

    // Plant state
    double y = 0.0;
    double setpoint = 1.0;
    constexpr double dt = 0.01;
    constexpr double a = 0.9;

    ctrlpp::implot::SignalRecorder recorder(4000);

    ctrlpp::implot::SimHarness sim(dt,
        [&](double t) {
            Vec sp = Vec::Constant(setpoint);
            Vec meas = Vec::Constant(y);
            auto u = pid.compute(sp, meas, dt);
            y = a * y + (1.0 - a) * u[0];

            recorder.record("setpoint", t, setpoint);
            recorder.record("measurement", t, y);
            recorder.record("control", t, u[0]);
        },
        [&] {
            pid.reset();
            y = 0.0;
            recorder.clear();
        });

    ctrlpp::implot::App app(1280, 720, "Basic PID Controller");

    app.run([&] {
        sim.advance(ImGui::GetIO().DeltaTime);

        ImGui::Begin("Simulation");
        sim.draw_controls();
        ImGui::End();

        ImGui::Begin("PID Parameters");

        auto sp_f = static_cast<float>(setpoint);
        if (ImGui::SliderFloat("Setpoint", &sp_f, -5.0f, 5.0f))
            setpoint = sp_f;

        ImGui::Separator();

        auto kp_f = static_cast<float>(cfg.kp[0]);
        auto ki_f = static_cast<float>(cfg.ki[0]);
        auto kd_f = static_cast<float>(cfg.kd[0]);
        bool changed = false;
        changed |= ImGui::SliderFloat("Kp", &kp_f, 0.0f, 20.0f);
        changed |= ImGui::SliderFloat("Ki", &ki_f, 0.0f, 10.0f);
        changed |= ImGui::SliderFloat("Kd", &kd_f, 0.0f, 5.0f);
        if (changed) {
            cfg.kp[0] = kp_f;
            cfg.ki[0] = ki_f;
            cfg.kd[0] = kd_f;
            pid.set_params(cfg);
        }

        ImGui::Separator();

        auto out_min_f = static_cast<float>(cfg.output_min[0]);
        auto out_max_f = static_cast<float>(cfg.output_max[0]);
        bool limits_changed = false;
        limits_changed |= ImGui::SliderFloat("Output Min", &out_min_f, -50.0f, 0.0f);
        limits_changed |= ImGui::SliderFloat("Output Max", &out_max_f, 0.0f, 50.0f);
        if (limits_changed) {
            cfg.output_min[0] = out_min_f;
            cfg.output_max[0] = out_max_f;
            pid.set_params(cfg);
        }

        if (ImGui::Button("Export SVG"))
            ctrlpp::implot::write_svg("basic_pid.svg", recorder,
                {"setpoint", "measurement", "control"});

        ImGui::End();

        auto t = static_cast<float>(sim.sim_time());
        ctrlpp::implot::scrolling_plot(
            {"Response", "Time (s)", "Value", 10.0f},
            recorder, {"setpoint", "measurement"}, t);
        ctrlpp::implot::scrolling_plot(
            {"Control Signal", "Time (s)", "u", 10.0f},
            recorder, {"control"}, t);
    });
}
