#include <ctrlpp/pid.h>
#include <ctrlpp/pid_policies.h>
#include <ctrlpp/eigen_linalg.h>

#include <cstdio>

int main()
{
    using Policy = ctrlpp::EigenLinalgPolicy;
    using Vec = Eigen::Matrix<double, 1, 1>;

    // Outer loop: temperature PI with anti-windup and IAE performance metric
    using OuterPid = ctrlpp::Pid<double, 1, 1, 1, Policy,
        ctrlpp::AntiWindup<ctrlpp::BackCalc>,
        ctrlpp::PerfAssessment<ctrlpp::IAE>>;

    OuterPid::config_type outer_cfg{};
    outer_cfg.kp = Vec::Constant(1.5);
    outer_cfg.ki = Vec::Constant(0.3);
    outer_cfg.kd = Vec::Constant(0.0);
    outer_cfg.output_min = Vec::Constant(0.0);
    outer_cfg.output_max = Vec::Constant(100.0);

    OuterPid outer(outer_cfg);

    // Inner loop: heater power PID with derivative filter
    using InnerPid = ctrlpp::Pid<double, 1, 1, 1, Policy,
        ctrlpp::DerivFilter>;

    InnerPid::config_type inner_cfg{};
    inner_cfg.kp = Vec::Constant(3.0);
    inner_cfg.ki = Vec::Constant(1.0);
    inner_cfg.kd = Vec::Constant(0.2);
    inner_cfg.output_min = Vec::Constant(0.0);
    inner_cfg.output_max = Vec::Constant(100.0);
    inner_cfg.template policy<ctrlpp::DerivFilter>().n = {10.0};

    InnerPid inner(inner_cfg);

    // Plant states
    double temperature = 20.0;
    double heater_power = 0.0;
    constexpr double temp_setpoint = 60.0;
    constexpr double dt = 0.1;
    constexpr int steps = 500;

    std::printf("Cascade PID temperature control\n");
    std::printf("Setpoint: %.1f deg (gain scheduling at step 200)\n\n", temp_setpoint);
    std::printf("%-6s  %-6s  %-8s  %-8s  %-10s  %-10s\n",
                "step", "time", "temp_sp", "temp", "heater_sp", "heater");
    std::printf("------  ------  --------  --------  ----------  ----------\n");

    for (int k = 0; k < steps; ++k) {
        // Gain scheduling: at step 200, increase outer gains (operating point change)
        if (k == 200) {
            outer_cfg.kp = Vec::Constant(2.0);
            outer_cfg.ki = Vec::Constant(0.5);
            outer.set_params(outer_cfg);
            std::printf(">>> Gain schedule change at step %d\n", k);
        }

        Vec sp_outer = Vec::Constant(temp_setpoint);
        Vec meas_outer = Vec::Constant(temperature);

        // Outer loop: temperature -> heater power setpoint
        // Use tracking form: outer tracks inner output for bumpless cascade
        Vec inner_output = Vec::Constant(heater_power);
        auto heater_sp_vec = outer.compute(sp_outer, meas_outer, dt, inner_output);
        double heater_sp = heater_sp_vec[0];

        // Inner loop: heater power setpoint -> control voltage
        Vec sp_inner = Vec::Constant(heater_sp);
        Vec meas_inner = Vec::Constant(heater_power);
        auto control = inner.compute(sp_inner, meas_inner, dt);

        // Print every 10th step
        if (k % 10 == 0) {
            std::printf("%-6d  %6.1f  %8.1f  %8.3f  %10.3f  %10.3f\n",
                        k, static_cast<double>(k) * dt, temp_setpoint,
                        temperature, heater_sp, heater_power);
        }

        // Plant dynamics
        // Inner: heater_power(k+1) = 0.8 * heater_power(k) + 0.2 * u_inner(k)
        heater_power = 0.8 * heater_power + 0.2 * control[0];
        // Outer: temperature(k+1) = 0.99 * temperature(k) + 0.01 * heater_power(k)
        temperature = 0.99 * temperature + 0.01 * heater_power;
    }

    // Print final IAE metric
    auto iae = outer.metric<ctrlpp::IAE>();
    std::printf("\n--- Performance ---\n");
    std::printf("Final temperature: %.3f (setpoint: %.1f)\n", temperature, temp_setpoint);
    std::printf("IAE metric:        %.4f\n", iae[0]);
}
