#include "ctrlpp/control/lqr.h"
#include "ctrlpp/control/dare.h"


#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

TEST_CASE("lqr_gain scalar integrator")
{
    // A=1, B=1, Q=1, R=1
    // DARE gives P = golden ratio, K = P / (1 + P)
    Eigen::Matrix<double, 1, 1> A, B, Q, R;
    A(0, 0) = 1.0;
    B(0, 0) = 1.0;
    Q(0, 0) = 1.0;
    R(0, 0) = 1.0;

    auto result = ctrlpp::lqr_gain<double, 1, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    double golden = (1.0 + std::sqrt(5.0)) / 2.0;
    double expected_K = golden / (1.0 + golden);
    CHECK_THAT((*result)(0, 0), Catch::Matchers::WithinAbs(expected_K, 1e-10));
}

TEST_CASE("lqr_gain double integrator stabilizes system")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, 1.0, 0.0, 1.0;
    B << 0.5, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;

    auto result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    auto K = *result;
    Eigen::Matrix<double, 2, 2> Acl = A - B * K;
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> solver(Acl, false);
    for(int i = 0; i < 2; ++i)
        CHECK(std::abs(solver.eigenvalues()(i)) < 1.0);
}

TEST_CASE("lqr_gain with cross-weight N")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B, N;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, 1.0, 0.0, 1.0;
    B << 0.5, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;
    N << 0.1, 0.2;

    auto result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R, N);
    REQUIRE(result.has_value());

    // Verify stabilizing
    auto K = *result;
    Eigen::Matrix<double, 2, 2> Acl = A - B * K;
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> solver(Acl, false);
    for(int i = 0; i < 2; ++i)
        CHECK(std::abs(solver.eigenvalues()(i)) < 1.0);
}

TEST_CASE("lqr_gain non-stabilizable returns nullopt")
{
    // A has unstable mode at eigenvalue 2, B cannot reach it
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 2.0, 0.0, 0.0, 0.5;
    B << 0.0, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;

    auto result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    CHECK_FALSE(result.has_value());
}

TEST_CASE("lqr_finite converges to infinite-horizon gain")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, 1.0, 0.0, 1.0;
    B << 0.5, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;

    auto Qf = Q;

    // Infinite-horizon gain for reference
    auto K_inf_opt = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    REQUIRE(K_inf_opt.has_value());
    auto K_inf = *K_inf_opt;

    // Long horizon: first gain should be close to infinite-horizon
    auto gains = ctrlpp::lqr_finite<double, 2, 1>(A, B, Q, R, Qf, 200);
    REQUIRE(gains.size() == 200);

    // K_0 should be close to K_inf for large horizon
    CHECK((gains[0] - K_inf).norm() < 1e-6);
}

TEST_CASE("lqr_tv_gains with constant matrices matches lqr_finite")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, 1.0, 0.0, 1.0;
    B << 0.5, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;
    auto Qf = Q;

    constexpr std::size_t horizon = 20;

    auto gains_finite = ctrlpp::lqr_finite<double, 2, 1>(A, B, Q, R, Qf, horizon);

    // Create constant vectors for time-varying interface
    std::vector<Eigen::Matrix<double, 2, 2>> As(horizon, A);
    std::vector<Eigen::Matrix<double, 2, 1>> Bs(horizon, B);
    std::vector<Eigen::Matrix<double, 2, 2>> Qs(horizon, Q);
    std::vector<Eigen::Matrix<double, 1, 1>> Rs(horizon, R);

    auto gains_tv = ctrlpp::lqr_tv_gains<double, 2, 1>(As, Bs, Qs, Rs, Qf, horizon);

    REQUIRE(gains_tv.size() == gains_finite.size());
    for(std::size_t k = 0; k < horizon; ++k)
        CHECK((gains_tv[k] - gains_finite[k]).norm() < 1e-12);
}

TEST_CASE("lqi_result achieves zero steady-state error")
{
    // First-order system: A=0.9, B=1, C=1
    // With integral action, should eliminate steady-state error to step reference
    Eigen::Matrix<double, 1, 1> A, B, C;
    A(0, 0) = 0.9;
    B(0, 0) = 1.0;
    C(0, 0) = 1.0;

    // Augmented Q and R
    Eigen::Matrix<double, 2, 2> Q_aug = Eigen::Matrix<double, 2, 2>::Identity();
    Q_aug(1, 1) = 10.0; // Higher weight on integral state
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 1.0;

    auto result = ctrlpp::lqi_gain<double, 1, 1, 1>(A, B, C, Q_aug, R);
    REQUIRE(result.has_value());

    auto Kx = result->Kx;
    auto Ki = result->Ki;

    // Verify closed-loop stability of augmented system
    // A_aug = [[A, 0], [-C, I]], B_aug = [[B], [0]]
    // K_aug = [Kx, Ki]
    // A_cl = A_aug - B_aug * K_aug
    Eigen::Matrix<double, 2, 2> A_aug;
    A_aug << A(0, 0), 0.0, -C(0, 0), 1.0;
    Eigen::Matrix<double, 2, 1> B_aug;
    B_aug << B(0, 0), 0.0;
    Eigen::Matrix<double, 1, 2> K_aug;
    K_aug << Kx(0, 0), Ki(0, 0);
    Eigen::Matrix<double, 2, 2> A_cl = A_aug - B_aug * K_aug;
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> eigsolver(A_cl, false);
    for(int i = 0; i < 2; ++i)
        REQUIRE(std::abs(eigsolver.eigenvalues()(i)) < 1.0);

    // Simulate closed loop with step reference r=1.0
    // Augmented dynamics: xi_{k+1} = xi_k + (r_k - C*x_k)
    // Control: u = -(Kx*x + Ki*xi)
    double x = 0.0;
    double xi = 0.0;
    double r = 1.0;

    for(int step = 0; step < 500; ++step)
    {
        double u = -(Kx(0, 0) * x + Ki(0, 0) * xi);
        double x_next = A(0, 0) * x + B(0, 0) * u;
        double y = C(0, 0) * x;
        xi = xi + (r - y);
        x = x_next;
    }

    double y_final = C(0, 0) * x;
    CHECK_THAT(y_final, Catch::Matchers::WithinAbs(r, 1e-4));
}

TEST_CASE("lqr_cost evaluates trajectory cost")
{
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 1.0;

    // Simple 2-step trajectory
    std::vector<Eigen::Matrix<double, 2, 1>> xs(3);
    std::vector<Eigen::Matrix<double, 1, 1>> us(2);

    xs[0] << 1.0, 0.0;
    xs[1] << 0.5, 0.1;
    xs[2] << 0.1, 0.05;

    us[0](0, 0) = 0.5;
    us[1](0, 0) = 0.3;

    // Expected: sum of x^T Q x + u^T R u for steps 0,1, plus x_2^T Q x_2 terminal
    double expected = 0.0;
    // Step 0: x0^T Q x0 + u0^T R u0 = 1.0 + 0.25 = 1.25
    expected += 1.0 * 1.0 + 0.0 * 0.0 + 0.5 * 0.5;
    // Step 1: x1^T Q x1 + u1^T R u1 = 0.25 + 0.01 + 0.09 = 0.35
    expected += 0.5 * 0.5 + 0.1 * 0.1 + 0.3 * 0.3;
    // Terminal: x2^T Q x2 = 0.01 + 0.0025 = 0.0125
    expected += 0.1 * 0.1 + 0.05 * 0.05;

    auto cost = ctrlpp::lqr_cost<double, 2, 1>(std::span<const Eigen::Matrix<double, 2, 1>>{xs}, std::span<const Eigen::Matrix<double, 1, 1>>{us}, Q, R);

    CHECK_THAT(cost, Catch::Matchers::WithinAbs(expected, 1e-12));
}

TEST_CASE("lqr class compute returns -K*x")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, 1.0, 0.0, 1.0;
    B << 0.5, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;

    auto K_opt = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    REQUIRE(K_opt.has_value());
    auto K = *K_opt;

    ctrlpp::lqr<double, 2, 1> lqr(K);

    Eigen::Matrix<double, 2, 1> x;
    x << 1.0, 0.5;

    auto u = lqr.compute(x);
    Eigen::Matrix<double, 1, 1> expected = -K * x;
    CHECK((u - expected).norm() < 1e-12);
    CHECK((lqr.gain() - K).norm() < 1e-12);
}

TEST_CASE("lqr_gain with N returns nullopt when DARE fails")
{
    // Non-stabilizable system with cross-weight N
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B, N;
    Eigen::Matrix<double, 1, 1> R;

    A << 2.0, 0.0, 0.0, 0.5;
    B << 0.0, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;
    N << 0.1, 0.2;

    auto result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R, N);
    CHECK_FALSE(result.has_value());
}

TEST_CASE("lqr_gain with singular A returns nullopt")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, 0.0, 0.0, 0.0; // singular
    B << 1.0, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;

    auto result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    CHECK_FALSE(result.has_value());
}

TEST_CASE("lqr_time_varying indexes correctly")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, 1.0, 0.0, 1.0;
    B << 0.5, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;
    auto Qf = Q;

    constexpr std::size_t horizon = 10;
    auto gains = ctrlpp::lqr_finite<double, 2, 1>(A, B, Q, R, Qf, horizon);

    ctrlpp::lqr_time_varying<double, 2, 1> lqr_tv(std::move(gains));
    CHECK(lqr_tv.horizon() == horizon);

    Eigen::Matrix<double, 2, 1> x;
    x << 1.0, -0.5;

    // Recompute gains for verification
    auto gains_ref = ctrlpp::lqr_finite<double, 2, 1>(A, B, Q, R, Qf, horizon);

    for(std::size_t k = 0; k < horizon; ++k)
    {
        auto u = lqr_tv.compute(x, k);
        Eigen::Matrix<double, 1, 1> expected = -gains_ref[k] * x;
        CHECK((u - expected).norm() < 1e-12);
        CHECK((lqr_tv.gain(k) - gains_ref[k]).norm() < 1e-12);
    }
}
