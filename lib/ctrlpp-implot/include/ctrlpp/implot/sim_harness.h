#ifndef HPP_GUARD_CTRLPP_IMPLOT_SIM_HARNESS_H
#define HPP_GUARD_CTRLPP_IMPLOT_SIM_HARNESS_H

#include <cstddef>
#include <functional>

namespace ctrlpp::implot {

class SimHarness {
public:
    SimHarness(double dt, std::function<void(double)> step_fn,
               std::function<void()> reset_fn);

    void advance(double real_dt);
    void draw_controls();

    [[nodiscard]] double sim_time() const;
    [[nodiscard]] bool running() const;
    void set_speed(double speed);

private:
    std::function<void(double)> step_fn_;
    std::function<void()> reset_fn_;
    double dt_;
    double sim_time_{0.0};
    double accumulator_{0.0};
    double speed_{1.0};
    bool running_{false};
    bool single_step_{false};
    static constexpr std::size_t kMaxStepsPerFrame{10};
};

}

#endif
