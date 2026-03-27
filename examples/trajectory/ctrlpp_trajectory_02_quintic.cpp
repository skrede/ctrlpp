// Usage: ./ctrlpp_trajectory_02_quintic | gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '-' using 1:2 with lines title 'position', '' using 1:3 with lines title 'velocity', '' using 1:4 with lines title 'acceleration'"
// Redirect: ./ctrlpp_trajectory_02_quintic > output.csv

#include <ctrlpp/trajectory/quintic_trajectory.h>

#include <iomanip>
#include <iostream>

int main()
{
    // 1-DOF joint move: 0 -> 3.14 rad over 3 seconds
    // Zero velocity and acceleration at both ends
    using Vec = Eigen::Matrix<double, 1, 1>;

    Vec const q0 = Vec::Constant(0.0);
    Vec const q1 = Vec::Constant(3.14);
    Vec const v0 = Vec::Zero();
    Vec const v1 = Vec::Zero();
    Vec const a0 = Vec::Zero();
    Vec const a1 = Vec::Zero();
    constexpr double duration = 3.0;

    auto const trajectory = ctrlpp::make_quintic_trajectory(q0, q1, v0, v1, a0, a1, duration);

    constexpr double dt = 0.01;
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
