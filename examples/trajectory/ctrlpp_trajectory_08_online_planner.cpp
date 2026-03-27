// Usage: ./ctrlpp_trajectory_08_online_planner | gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '-' using 1:2 with lines title 'position', '' using 1:4 with lines title 'target'"
// Redirect: ./ctrlpp_trajectory_08_online_planner > output.csv

#include <ctrlpp/trajectory/online_planner_2nd.h>

#include <iomanip>
#include <iostream>

int main()
{
    // 2nd-order online planner tracking a sequence of target changes.
    // Simulates real-time joystick or sensor-driven commands.
    ctrlpp::online_planner_2nd<double> planner({.v_max = 2.0, .a_max = 5.0});

    constexpr double dt = 0.005;
    constexpr double total_time = 10.0;

    // Target change schedule
    double current_target = 5.0;
    planner.update(current_target);

    std::cout << "time,position,velocity,target\n";

    for (double t = 0.0; t <= total_time; t += dt)
    {
        // Change target at specified times
        if (t >= 6.0 && current_target != 8.0)
        {
            current_target = 8.0;
            planner.update(current_target);
        }
        else if (t >= 3.0 && t < 6.0 && current_target != 2.0)
        {
            current_target = 2.0;
            planner.update(current_target);
        }

        auto const pt = planner.sample(t);

        std::cout << std::fixed << std::setprecision(4) << t << "," << pt.position[0] << "," << pt.velocity[0] << "," << current_target << "\n";
    }
}
