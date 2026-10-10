#include "switch_shapes.h"
#include "switch_planners.h"
#include "switch_trapezoid.h"
#include "switch_third_order.h"
#include "envelope_legs.h"
#include "agreement_legs.h"
#include "acceleration_legs.h"
#include "precondition_legs.h"

#include "ctrlpp/trajectory/double_s_trajectory.h"
#include "ctrlpp/trajectory/trapezoidal_trajectory.h"

#include <cstdio>
#include <cstring>

namespace planner_overshoot_sweep
{

void run_limit_legs(bool print_points)
{
    std::printf("\n\nThird-order planner retargeted from its cruise, acceleration-ramp and acceleration-nulling\n");
    std::printf("windows: every sample of the new plan against v_max and a_max, and every dense pair of\n");
    std::printf("samples against j_max plus the counted rounding allowance.\n");
    run_retarget_leg<double>("Retarget envelope", false, print_points);
    run_retarget_leg<float>("Retarget envelope", false, print_points);
    run_retarget_leg<double>("Reversing", true, print_points);
    run_retarget_leg<float>("Reversing", true, print_points);

    std::printf("\n\nShape constructors handed a boundary speed stepped across v_max one representable value\n");
    std::printf("at a time, as v0 and as v1, in both frames, under each boundary-velocity rounding count.\n");
    run_precondition_leg<double_s_trajectory<double>>("double_s_trajectory", print_points);
    run_precondition_leg<double_s_trajectory<float>>("double_s_trajectory", print_points);
    run_precondition_leg<trapezoidal_trajectory<double>>("trapezoidal_trajectory", print_points);
    run_precondition_leg<trapezoidal_trajectory<float>>("trapezoidal_trajectory", print_points);
}

void run_switch_legs(bool print_points)
{
    std::printf("\n\nEvery representable instant within %d steps of each phase switch of the two shapes and the two\n",
                switch_window_steps);
    std::printf("online planners, over seeded boundary fractions, end velocities at the limit and profiles whose jerk\n");
    std::printf("phases are near or below the sample-time resolution: |v| against v_max, |a| against the bound.\n");
    run_trapezoid_switch_leg<double>(print_points);
    run_trapezoid_switch_leg<float>(print_points);
    run_double_s_switch_leg<double>(print_points);
    run_double_s_switch_leg<float>(print_points);
    run_second_order_switch_leg<double>(print_points);
    run_second_order_switch_leg<float>(print_points);
    run_third_order_switch_leg<double>(print_points);
    run_third_order_switch_leg<float>(print_points);
}

}

using planner_overshoot_sweep::run_limit_legs;
using planner_overshoot_sweep::run_switch_legs;
using planner_overshoot_sweep::run_agreement_legs;
using planner_overshoot_sweep::run_acceleration_legs;

auto main(int argc, char** argv) -> int
{
    bool print_points = true;
    for (int index = 1; index < argc; ++index) {
        if (std::strcmp(argv[index], "--summary") == 0) {
            print_points = false;
        }
    }

    std::printf("Planner overshoot verdict against the minimum displacement a carry-velocity"
                " shape needs\n");
    run_agreement_legs<double>(print_points);
    std::printf("\n\n");
    run_agreement_legs<float>(print_points);
    run_acceleration_legs(print_points);
    run_limit_legs(print_points);
    run_switch_legs(print_points);

    std::printf("\n");
    return 0;
}
