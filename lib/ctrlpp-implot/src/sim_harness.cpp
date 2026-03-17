#include "ctrlpp/implot/sim_harness.h"

#include "imgui.h"

#include <algorithm>

namespace ctrlpp::implot {

SimHarness::SimHarness(double dt, std::function<void(double)> step_fn,
                       std::function<void()> reset_fn)
    : step_fn_{std::move(step_fn)}
    , reset_fn_{std::move(reset_fn)}
    , dt_{dt}
{
}

void SimHarness::advance(double real_dt)
{
    if (!running_ && !single_step_) {
        return;
    }

    double scaled_dt = real_dt * speed_;
    if (single_step_) {
        scaled_dt = dt_;
        single_step_ = false;
    }

    accumulator_ += scaled_dt;
    accumulator_ =
        std::min(accumulator_, static_cast<double>(kMaxStepsPerFrame) * dt_);

    while (accumulator_ >= dt_) {
        step_fn_(dt_);
        sim_time_ += dt_;
        accumulator_ -= dt_;
    }
}

void SimHarness::draw_controls()
{
    if (ImGui::Button(running_ ? "Pause" : "Play")) {
        running_ = !running_;
    }
    ImGui::SameLine();
    if (ImGui::Button("Step")) {
        single_step_ = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("Reset")) {
        reset_fn_();
        sim_time_ = 0.0;
        accumulator_ = 0.0;
        running_ = false;
    }

    auto speed_f = static_cast<float>(speed_);
    ImGui::SliderFloat("Speed", &speed_f, 0.1f, 10.0f, "%.1fx");
    speed_ = static_cast<double>(speed_f);

    ImGui::Text("Time: %.3f s", sim_time_);
}

double SimHarness::sim_time() const
{
    return sim_time_;
}

bool SimHarness::running() const
{
    return running_;
}

void SimHarness::set_speed(double speed)
{
    speed_ = speed;
}

}
