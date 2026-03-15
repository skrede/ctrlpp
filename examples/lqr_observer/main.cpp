// LQR with Observer Example
//
// Demonstrates external composition of observers and LQR controller:
// - KalmanFilter and LuenbergerObserver on the same discrete double-integrator plant
// - Concept-based polymorphism via ObserverPolicy
// - observer.state() fed to lqr.compute() (external composition pattern)

#include "ctrlpp/ctrlpp.h"
#include "ctrlpp/ctrlpp_eigen.h"

#include <iomanip>
#include <iostream>
#include <string>

// Dimensions: 2 states (position, velocity), 1 input (force), 1 output (position)
constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr std::size_t NY = 1;

using Scalar = double;
using Mat2 = Eigen::Matrix<Scalar, 2, 2>;
using Vec2 = Eigen::Matrix<Scalar, 2, 1>;
using Mat1 = Eigen::Matrix<Scalar, 1, 1>;
using MatC = Eigen::Matrix<Scalar, 1, 2>;
using MatB = Eigen::Matrix<Scalar, 2, 1>;
using MatD = Eigen::Matrix<Scalar, 1, 1>;

using System = ctrlpp::DiscreteStateSpace<ctrlpp::EigenLinalgPolicy, Scalar, NX, NU, NY>;

// Generic simulation loop accepting any ObserverPolicy observer.
// Demonstrates concept-based interchangeability: the same function works with
// KalmanFilter and LuenbergerObserver without runtime polymorphism.
template<ctrlpp::ObserverPolicy Observer>
void run_simulation(
    const System& sys,
    const ctrlpp::Lqr<Scalar, NX, NU>& lqr,
    Observer& observer,
    const std::string& name,
    const Vec2& x0_true,
    std::size_t steps)
{
    std::cout << "\n=== " << name << " ===\n\n";
    std::cout << std::setw(5) << "step"
              << std::setw(12) << "x1_true"
              << std::setw(12) << "x2_true"
              << std::setw(12) << "x1_est"
              << std::setw(12) << "x2_est"
              << std::setw(12) << "u"
              << "\n";
    std::cout << std::string(65, '-') << "\n";

    Vec2 x_true = x0_true;
    Eigen::Matrix<Scalar, 1, 1> u_ctrl = Eigen::Matrix<Scalar, 1, 1>::Zero();

    for (std::size_t k = 0; k < steps; ++k) {
        // 1. Observer prediction step with previous control input
        observer.predict(u_ctrl);

        // 2. Simulate true plant dynamics
        x_true = (sys.A * x_true + sys.B * u_ctrl).eval();

        // 3. Measure output from true state
        Eigen::Matrix<Scalar, 1, 1> y = (sys.C * x_true).eval();

        // 4. Observer measurement update
        observer.update(y);

        // 5. External composition: observer estimate feeds LQR controller
        u_ctrl = lqr.compute(observer.state());

        // 6. Print state
        if (k % 5 == 0 || k < 5) {
            const auto& x_est = observer.state();
            std::cout << std::setw(5) << k
                      << std::fixed << std::setprecision(4)
                      << std::setw(12) << x_true(0)
                      << std::setw(12) << x_true(1)
                      << std::setw(12) << x_est(0)
                      << std::setw(12) << x_est(1)
                      << std::setw(12) << u_ctrl(0)
                      << "\n";
        }
    }

    const auto& x_est = observer.state();
    std::cout << "\nFinal true state:     [" << x_true(0) << ", " << x_true(1) << "]\n";
    std::cout << "Final estimated state: [" << x_est(0) << ", " << x_est(1) << "]\n";
    std::cout << "Estimation error norm: " << (x_true - x_est).norm() << "\n";
}

int main()
{
    // -----------------------------------------------------------------------
    // 1. Define discrete-time double integrator (dt = 0.1 s)
    //    x = [position, velocity]
    //    x_{k+1} = A x_k + B u_k
    //    y_k     = C x_k
    // -----------------------------------------------------------------------
    constexpr Scalar dt = 0.1;

    Mat2 A;
    A << 1.0, dt,
         0.0, 1.0;

    MatB B;
    B << 0.5 * dt * dt,
         dt;

    MatC C;
    C << 1.0, 0.0;

    MatD D = MatD::Zero();

    System sys{A, B, C, D};

    std::cout << "LQR with Observer Example\n";
    std::cout << "=========================\n\n";
    std::cout << "Plant: discrete double integrator (dt=" << dt << "s)\n";
    std::cout << "States: [position, velocity], Input: [force], Output: [position]\n\n";

    // -----------------------------------------------------------------------
    // 2. Check controllability and observability
    // -----------------------------------------------------------------------
    bool controllable = ctrlpp::is_controllable<ctrlpp::EigenLinalgPolicy, Scalar, NX, NU>(A, B);
    bool observable = ctrlpp::is_observable<ctrlpp::EigenLinalgPolicy, Scalar, NX, NY>(A, C);

    std::cout << "Controllability: " << (controllable ? "YES" : "NO") << "\n";
    std::cout << "Observability:   " << (observable ? "YES" : "NO") << "\n\n";

    // -----------------------------------------------------------------------
    // 3. Design LQR gain: Q = diag(10, 1), R = 1
    // -----------------------------------------------------------------------
    Mat2 Q = Mat2::Zero();
    Q(0, 0) = 10.0;
    Q(1, 1) = 1.0;

    Mat1 R;
    R << 1.0;

    auto K_opt = ctrlpp::lqr_gain<Scalar, NX, NU>(A, B, Q, R);
    if (!K_opt) {
        std::cerr << "LQR design failed (non-stabilizable system)\n";
        return 1;
    }

    ctrlpp::Lqr<Scalar, NX, NU> lqr(*K_opt);

    std::cout << "LQR gain K = [" << (*K_opt)(0, 0) << ", " << (*K_opt)(0, 1) << "]\n";

    // Verify closed-loop stability
    bool cl_stable = ctrlpp::is_stable_closed_loop<ctrlpp::EigenLinalgPolicy, Scalar, NX, NU>(
        A, B, *K_opt);
    std::cout << "Closed-loop stable: " << (cl_stable ? "YES" : "NO") << "\n";

    // -----------------------------------------------------------------------
    // 4. Initial conditions: true state starts at [1, 0], estimate at [0, 0]
    // -----------------------------------------------------------------------
    Vec2 x0_true;
    x0_true << 1.0, 0.0;

    Vec2 x0_est = Vec2::Zero();

    constexpr std::size_t sim_steps = 50;

    // -----------------------------------------------------------------------
    // 5. Run with Kalman filter
    //    Process noise Q_proc = 0.01*I, measurement noise R_meas = 0.1
    // -----------------------------------------------------------------------
    {
        Mat2 Q_proc = 0.01 * Mat2::Identity();
        Mat1 R_meas;
        R_meas << 0.1;
        Mat2 P0 = Mat2::Identity();

        ctrlpp::KalmanFilter<Scalar, NX, NU, NY> kf(sys, Q_proc, R_meas, x0_est, P0);

        run_simulation(sys, lqr, kf, "Kalman Filter Observer", x0_true, sim_steps);
    }

    // -----------------------------------------------------------------------
    // 6. Run with Luenberger observer
    //    Observer poles at 0.3 and 0.4 (faster than open-loop plant poles)
    // -----------------------------------------------------------------------
    {
        std::array<std::complex<Scalar>, NX> obs_poles = {
            std::complex<Scalar>{0.3, 0.0},
            std::complex<Scalar>{0.4, 0.0}
        };

        auto L_opt = ctrlpp::place_observer<Scalar, NX, NY>(A, C, obs_poles);
        if (!L_opt) {
            std::cerr << "Observer pole placement failed\n";
            return 1;
        }

        std::cout << "\nLuenberger observer gain L = [" << (*L_opt)(0, 0)
                  << ", " << (*L_opt)(1, 0) << "]^T\n";

        // Verify observer stability
        bool obs_stable = ctrlpp::is_stable_observer<ctrlpp::EigenLinalgPolicy, Scalar, NX, NY>(
            A, *L_opt, C);
        std::cout << "Observer stable: " << (obs_stable ? "YES" : "NO") << "\n";

        ctrlpp::LuenbergerObserver<Scalar, NX, NU, NY> luenberger(sys, *L_opt, x0_est);

        run_simulation(sys, lqr, luenberger, "Luenberger Observer", x0_true, sim_steps);
    }

    std::cout << "\nBoth observers converge to the true state while LQR regulates to zero.\n";
    std::cout << "The same run_simulation() template works with either observer via ObserverPolicy.\n";

    return 0;
}
