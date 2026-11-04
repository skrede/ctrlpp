// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_trajectory_08_online_planner' using 1:2 with lines title 'position', '' using 1:4 with lines title 'target'"
// Redirect: ./ctrlpp_trajectory_08_online_planner > output.csv

#include "ctrlpp/trajectory/online_planner_2nd.h"

#include <iomanip>
#include <iostream>

int main()
{
    // The retargets stand in for joystick or sensor commands arriving in real time. create refuses limits
    // that are not finite and strictly positive, since the planner divides by them.
    auto planner_result = ctrlpp::online_planner_2nd<double>::create({.v_max = 2.0, .a_max = 5.0});
    if(!planner_result.has_value())
    {
        std::cerr << "invalid planner limits\n";
        return 1;
    }
    auto& planner = *planner_result;

    constexpr double dt = 0.005;
    constexpr double total_time = 10.0;

    double current_target = 5.0;
    if(!planner.update(current_target).has_value())
    {
        std::cerr << "invalid target\n";
        return 1;
    }

    // A reversal or overshoot mid-motion is braked to rest and replanned from the stopping point, which
    // diagnostics() reports; the limits hold either way, only the move takes longer.
    std::cout << "time,position,velocity,target,substituted\n";

    for(double t = 0.0; t <= total_time; t += dt)
    {
        if(t >= 6.0 && current_target != 8.0)
        {
            current_target = 8.0;
            if(!planner.update(current_target).has_value())
                return 1;
        }
        else if(t >= 3.0 && t < 6.0 && current_target != 2.0)
        {
            current_target = 2.0;
            if(!planner.update(current_target).has_value())
                return 1;
        }

        auto const pt = planner.sample(t);

        bool const substituted = planner.diagnostics().disposition == ctrlpp::online_planner_disposition::braked_and_replanned;

        std::cout << std::fixed << std::setprecision(4) << t << "," << pt.position[0] << "," << pt.velocity[0] << "," << current_target << "," << (substituted ? 1 : 0) << "\n";
    }
}
