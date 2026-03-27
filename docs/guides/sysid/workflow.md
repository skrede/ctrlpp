# System Identification Workflow

System identification turns measured input-output data into a mathematical
model. ctrlpp provides the full pipeline: collect data, identify a model,
and use that model for control design.

## The Three Phases

1. **Collect data** -- apply excitation signals and record the response
2. **Identify model** -- fit an ARX or state-space model to the data
3. **Use model** -- feed the identified `discrete_state_space` into LQR, MPC,
   or a Kalman filter

## Phase 1: Collect Data

Generate a pseudo-random excitation signal and record the plant response. In
practice this data comes from real sensors; here we simulate a first-order
system:

```
y(t) = 0.7 * y(t-1) + 0.3 * u(t-1)
```

## Phase 2: Identify the Model

### Batch ARX

Batch ARX solves a least-squares problem over the full dataset. It returns a
`discrete_state_space` in observer canonical form plus fit metrics.

### Online RLS

Recursive least squares (RLS) identifies parameters online, one sample at a
time. It adapts to slowly varying systems via the forgetting factor.

## Complete Program: Batch ARX to LQR

This example identifies a system from data and immediately uses the result
for LQR control design.

```cpp
// Usage: ./sysid_workflow | gnuplot -p -e "set datafile separator ','; plot '-' skip 1 using 1:2 with lines title 'state', '' using 1:3 with lines title 'control'"
#include <ctrlpp/sysid.h>
#include <ctrlpp/control/lqr.h>
#include <ctrlpp/model/propagate.h>

#include <Eigen/Dense>

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <random>

int main()
{
    // --- Phase 1: Generate data ---
    constexpr std::size_t N = 500;

    Eigen::Matrix<double, 1, static_cast<int>(N)> Y;
    Eigen::Matrix<double, 1, static_cast<int>(N)> U;

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> u_dist(-1.0, 1.0);

    double y = 0.0;
    double u_prev = 0.0;
    for (std::size_t t = 0; t < N; ++t)
    {
        double u = u_dist(gen);
        double y_new = 0.7 * y + 0.3 * u_prev;
        Y(0, static_cast<int>(t)) = y_new;
        U(0, static_cast<int>(t)) = u;
        y = y_new;
        u_prev = u;
    }

    // --- Phase 2: Identify model ---
    auto result = ctrlpp::batch_arx<1, 1>(Y, U);

    std::cerr << "NRMSE=" << result.metrics.nrmse
              << " VAF=" << result.metrics.vaf << "%\n";

    // --- Phase 3: Use model for control ---
    auto const& sys = result.system;

    // Design LQR for the identified system
    constexpr std::size_t NX = 1;
    constexpr std::size_t NU_ctrl = 1;

    Eigen::Matrix<double, NX, NX> Q_lqr = Eigen::Matrix<double, 1, 1>::Identity() * 10.0;
    Eigen::Matrix<double, NU_ctrl, NU_ctrl> R_lqr = Eigen::Matrix<double, 1, 1>::Identity();

    auto K_opt = ctrlpp::lqr_gain<double, NX, NU_ctrl>(sys.A, sys.B, Q_lqr, R_lqr);

    if (!K_opt)
    {
        std::cerr << "LQR design failed\n";
        return 1;
    }

    ctrlpp::lqr<double, NX, NU_ctrl> controller(*K_opt);

    // Closed-loop simulation with identified model
    Eigen::Vector<double, NX> x = Eigen::Vector<double, NX>::Constant(1.0);

    std::cout << "step,state,control\n";
    for (int k = 0; k < 50; ++k)
    {
        auto u_ctrl = controller.compute(x);
        std::cout << k << "," << x(0) << "," << u_ctrl(0) << "\n";
        x = sys.A * x + sys.B * u_ctrl;
    }
}
```

## Complete Program: Online RLS

For systems with slowly varying parameters, use RLS to track changes online:

```cpp
// Usage: ./rls_online | gnuplot -p -e "set datafile separator ','; plot '-' skip 1 using 1:2 with lines title 'param 1', '' using 1:3 with lines title 'param 2'"
#include <ctrlpp/sysid.h>

#include <Eigen/Dense>

#include <iostream>
#include <random>

int main()
{
    constexpr std::size_t NP = 2;
    constexpr double lambda = 0.98;

    ctrlpp::rls_config<double, NP> cfg;
    cfg.lambda = lambda;
    ctrlpp::rls<double, NP> estimator(cfg);

    Eigen::Vector2d true_params;
    true_params << 2.0, 0.5;

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> input(-1.0, 1.0);
    std::normal_distribution<double> noise(0.0, 0.01);

    std::cout << "step,param_1,param_2\n";

    for (int t = 1; t <= 200; ++t)
    {
        Eigen::Vector2d phi;
        phi << input(gen), input(gen);
        double y_meas = true_params.dot(phi) + noise(gen);

        estimator.update(y_meas, phi);

        auto theta = estimator.parameters();
        std::cout << t << "," << theta(0) << "," << theta(1) << "\n";
    }
}
```

## Choosing an Identification Method

| Method       | Data        | Use case                              |
| ------------ | ----------- | ------------------------------------- |
| `batch_arx`  | Offline     | Known model order, full dataset       |
| `n4sid`       | Offline     | Unknown order, subspace method        |
| `rls`         | Online      | Scalar output, time-varying params    |
| `recursive_arx`| Online    | ARX structure, streaming data         |

All offline methods return a `discrete_state_space` that plugs directly into
LQR, MPC, or Kalman filter constructors.

## Next Steps

- [Batch ARX API](../../api/sysid/batch-arx.md) -- offline identification
- [N4SID API](../../api/sysid/n4sid.md) -- subspace identification
- [RLS API](../../api/sysid/rls.md) -- recursive least squares
- [Fit Metrics API](../../api/sysid/fit-metrics.md) -- NRMSE, VAF
- [State Space API](../../api/model/state-space.md) -- model representation
- [Sysid Theory](../../background/sysid.md) -- ARX, subspace,
  recursive estimation theory
