// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_trajectory_03_trapezoidal' using 1:2 with lines title 'position', '' using 1:3 with lines title 'velocity', '' using 1:4 with lines title 'acceleration'"
// Redirect: ./ctrlpp_trajectory_03_trapezoidal > output.csv

#include <ctrlpp/trajectory/trapezoidal_trajectory.h>

#include <iomanip>
#include <iostream>

int main()
{
    // CNC linear axis move: 0 -> 100 mm
    // v_max = 50 mm/s, a_max = 200 mm/s^2, starting and ending at rest
    //
    // create() is the only way to build a profile, and it hands back the reason
    // when the command has none. A profile object therefore always satisfies its
    // own contract, so evaluate() below never has to be second-guessed.
    auto const created = ctrlpp::trapezoidal_trajectory<double>::create({
        .q0 = 0.0,
        .q1 = 100.0,
        .v_max = 50.0,
        .a_max = 200.0,
    });
    if (!created)
    {
        std::cerr << "The commanded move has no trapezoidal profile\n";
        return 1;
    }
    auto const& trajectory = created.value();

    constexpr double dt = 0.001;
    std::cout << "time,position,velocity,acceleration\n";

    for (double t = 0.0; t <= trajectory.duration(); t += dt)
    {
        auto const pt = trajectory.evaluate(t);
        std::cout << std::fixed << std::setprecision(4)
                  << t << ","
                  << pt.position[0] << "," << pt.velocity[0] << ","
                  << pt.acceleration[0] << "\n";
    }
}
