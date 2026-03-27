// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_trajectory_01_cubic' using 1:2 with lines title 'joint1', '' using 1:4 with lines title 'joint2'"
// Redirect: ./ctrlpp_trajectory_01_cubic > output.csv

#include <ctrlpp/trajectory/cubic_trajectory.h>

#include <iomanip>
#include <iostream>

int main()
{
    // 2-DOF robot arm joint move: (0, 0) -> (1.57, 0.78) rad over 2 seconds
    // Starting and ending at rest (v0 = v1 = 0)
    Eigen::Vector2d const q0{0.0, 0.0};
    Eigen::Vector2d const q1{1.57, 0.78};
    Eigen::Vector2d const v0{0.0, 0.0};
    Eigen::Vector2d const v1{0.0, 0.0};
    constexpr double duration = 2.0;

    auto const trajectory = ctrlpp::make_cubic_trajectory(q0, q1, v0, v1, duration);

    constexpr double dt = 0.01;
    std::cout << "time,q1_pos,q1_vel,q2_pos,q2_vel\n";

    for (double t = 0.0; t <= trajectory.duration(); t += dt)
    {
        auto const pt = trajectory.evaluate(t);
        std::cout << std::fixed << std::setprecision(4)
                  << t << ","
                  << pt.position[0] << "," << pt.velocity[0] << ","
                  << pt.position[1] << "," << pt.velocity[1] << "\n";
    }
}
