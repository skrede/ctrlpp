// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_trajectory_04_double_s' using 1:2 with lines title 'position', '' using 1:3 with lines title 'velocity', '' using 1:4 with lines title 'acceleration'"
// Redirect: ./ctrlpp_trajectory_04_double_s > output.csv

#include "ctrlpp/trajectory/double_s_trajectory.h"

#include <iomanip>
#include <iostream>

int main()
{
    // A joint in rad, rad/s, rad/s^2 and rad/s^3.
    // create() is the only way to build a profile and reports why a command has none, so a profile always
    // satisfies its own contract and evaluate() needs no second-guessing.
    auto const created = ctrlpp::double_s_trajectory<double>::create({
            .q0 = 0.0,
            .q1 = 2.0,
            .v_max = 1.0,
            .a_max = 5.0,
            .j_max = 50.0,
    });
    if(!created)
    {
        std::cerr << "The commanded move has no double-S profile\n";
        return 1;
    }
    auto const& trajectory = *created;

    constexpr double dt = 0.001;
    std::cout << "time,position,velocity,acceleration\n";

    for(double t = 0.0; t <= trajectory.duration(); t += dt)
    {
        auto const pt = trajectory.evaluate(t);
        std::cout << std::fixed << std::setprecision(4) << t << "," << pt.position[0] << "," << pt.velocity[0] << "," << pt.acceleration[0] << "\n";
    }
}
