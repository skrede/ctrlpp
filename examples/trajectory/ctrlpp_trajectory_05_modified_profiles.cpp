// Usage: ./ctrlpp_trajectory_05_modified_profiles | gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '-' using 1:4 with lines title 'trap_acc', '' using 1:7 with lines title 'sin_acc'"
// Redirect: ./ctrlpp_trajectory_05_modified_profiles > output.csv

#include <ctrlpp/trajectory/modified_sin_trajectory.h>
#include <ctrlpp/trajectory/modified_trap_trajectory.h>

#include <iomanip>
#include <iostream>

int main()
{
    // Vibration-sensitive cam design: 10 mm move in 0.5 seconds
    // Compare modified trapezoidal vs modified sinusoidal profiles
    ctrlpp::modified_trap_trajectory<double> trap({.q0 = 0.0, .q1 = 10.0, .T = 0.5});
    ctrlpp::modified_sin_trajectory<double> sin_prof({.q0 = 0.0, .q1 = 10.0, .T = 0.5});

    constexpr double dt = 0.001;
    constexpr double T = 0.5;
    std::cout << "time,trap_pos,trap_vel,trap_acc,sin_pos,sin_vel,sin_acc\n";

    for (double t = 0.0; t <= T; t += dt)
    {
        auto const tp = trap.evaluate(t);
        auto const sp = sin_prof.evaluate(t);
        std::cout << std::fixed << std::setprecision(4)
                  << t << ","
                  << tp.position[0] << "," << tp.velocity[0] << "," << tp.acceleration[0] << ","
                  << sp.position[0] << "," << sp.velocity[0] << "," << sp.acceleration[0] << "\n";
    }
}
