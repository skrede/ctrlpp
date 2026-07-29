// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_trajectory_08_online_planner' using 1:2 with lines title 'position', '' using 1:4 with lines title 'target'"
// Redirect: ./ctrlpp_trajectory_08_online_planner > output.csv

#include <ctrlpp/trajectory/online_planner_2nd.h>

#include <iomanip>
#include <iostream>

int main()
{
    // 2nd-order online planner tracking a sequence of target changes.
    // Simulates real-time joystick or sensor-driven commands.
    // create validates the kinematic limits (finite and strictly positive,
    // because they divide in the planner math) and reports rejections through
    // ctrlpp::expected; unwrap after checking.
    auto planner_result =
        ctrlpp::online_planner_2nd<double>::create({.v_max = 2.0, .a_max = 5.0});
    if (!planner_result.has_value())
    {
        std::cerr << "invalid planner limits\n";
        return 1;
    }
    auto& planner = *planner_result;

    constexpr double dt = 0.005;
    constexpr double total_time = 10.0;

    // Target change schedule
    double current_target = 5.0;
    if (!planner.update(current_target).has_value())
    {
        std::cerr << "invalid target\n";
        return 1;
    }

    // A retarget mid-motion does not always produce the commanded profile: a
    // reversal or an overshoot is braked to rest and replanned from the stopping
    // point instead. A finite accepted target still reports the resulting
    // disposition through diagnostics(). The motion respects the same limits
    // either way; what changes is how long the move takes.
    std::cout << "time,position,velocity,target,substituted\n";

    for (double t = 0.0; t <= total_time; t += dt)
    {
        // Change target at specified times
        if (t >= 6.0 && current_target != 8.0)
        {
            current_target = 8.0;
            if (!planner.update(current_target).has_value())
                return 1;
        }
        else if (t >= 3.0 && t < 6.0 && current_target != 2.0)
        {
            current_target = 2.0;
            if (!planner.update(current_target).has_value())
                return 1;
        }

        auto const pt = planner.sample(t);

        bool const substituted = planner.diagnostics().disposition == ctrlpp::online_planner_disposition::braked_and_replanned;

        std::cout << std::fixed << std::setprecision(4) << t << "," << pt.position[0] << "," << pt.velocity[0] << "," << current_target << "," << (substituted ? 1 : 0) << "\n";
    }
}
