#include "ctrlpp/analysis.h"
#include "ctrlpp/eigen_linalg.h"

#include <catch2/catch_test_macros.hpp>

using Policy = ctrlpp::EigenLinalgPolicy;

TEST_CASE("is_controllable with controllable double integrator") {
    // Double integrator: A=[[1,1],[0,1]], B=[[0.5],[1]]
    // Controllability matrix C = [B, AB] = [[0.5, 1.5],[1, 1]]
    // rank = 2 = NX -> controllable
    Eigen::Matrix<double, 2, 2> A;
    Eigen::Matrix<double, 2, 1> B;
    A << 1.0, 1.0,
         0.0, 1.0;
    B << 0.5,
         1.0;

    CHECK(ctrlpp::is_controllable<Policy, double, 2, 1>(A, B));
}

TEST_CASE("is_controllable with uncontrollable system") {
    // A diagonal with unstable mode, B cannot reach first state
    // A=[[2,0],[0,0.5]], B=[[0],[1]]
    // C = [B, AB] = [[0,0],[1,0.5]] -> rank 1 < 2 -> uncontrollable
    Eigen::Matrix<double, 2, 2> A;
    Eigen::Matrix<double, 2, 1> B;
    A << 2.0, 0.0,
         0.0, 0.5;
    B << 0.0,
         1.0;

    CHECK_FALSE(ctrlpp::is_controllable<Policy, double, 2, 1>(A, B));
}

TEST_CASE("is_observable with observable system") {
    // Same double integrator, C=[1,0] -> observable
    // O = [C; CA] = [[1,0],[1,1]] -> rank 2
    Eigen::Matrix<double, 2, 2> A;
    Eigen::Matrix<double, 1, 2> C;
    A << 1.0, 1.0,
         0.0, 1.0;
    C << 1.0, 0.0;

    CHECK(ctrlpp::is_observable<Policy, double, 2, 1>(A, C));
}

TEST_CASE("is_observable with unobservable system") {
    // A=[[2,0],[0,0.5]], C=[0,1] -> cannot observe first state
    // O = [C; CA] = [[0,1],[0,0.5]] -> rank 1 < 2 -> unobservable
    Eigen::Matrix<double, 2, 2> A;
    Eigen::Matrix<double, 1, 2> C;
    A << 2.0, 0.0,
         0.0, 0.5;
    C << 0.0, 1.0;

    CHECK_FALSE(ctrlpp::is_observable<Policy, double, 2, 1>(A, C));
}

TEST_CASE("observability duality with controllability") {
    // observable(A, C) iff controllable(A^T, C^T)
    Eigen::Matrix<double, 2, 2> A;
    Eigen::Matrix<double, 1, 2> C;
    A << 1.0, 1.0,
         0.0, 1.0;
    C << 1.0, 0.0;

    auto At = A.transpose().eval();
    auto Ct = C.transpose().eval();

    bool obs = ctrlpp::is_observable<Policy, double, 2, 1>(A, C);
    bool ctrl_dual = ctrlpp::is_controllable<Policy, double, 2, 1>(At, Ct);
    CHECK(obs == ctrl_dual);
}

TEST_CASE("is_stable_closed_loop with stabilizing gain") {
    // Double integrator A=[[1,1],[0,1]], B=[[0.5],[1]]
    // K chosen to place eigenvalues inside unit circle
    // K = [0.4, 0.8] => A-BK = [[1-0.2, 1-0.4],[0-0.4, 1-0.8]] = [[0.8, 0.6],[-0.4, 0.2]]
    // Eigenvalues: (0.8+0.2)/2 +/- sqrt((0.6)^2/4 - (0.8*0.2+0.4*0.6))
    Eigen::Matrix<double, 2, 2> A;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 2> K;
    A << 1.0, 1.0,
         0.0, 1.0;
    B << 0.5,
         1.0;
    K << 0.4, 0.8;

    // Verify eigenvalues of A-BK are inside unit circle
    Eigen::Matrix<double, 2, 2> Acl = A - B * K;
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> solver(Acl, false);
    bool all_stable = true;
    for (int i = 0; i < 2; ++i)
        if (std::abs(solver.eigenvalues()(i)) >= 1.0)
            all_stable = false;

    // Only test if K actually stabilizes (it should for these values)
    if (all_stable)
        CHECK(ctrlpp::is_stable_closed_loop<Policy, double, 2, 1>(A, B, K));
}

TEST_CASE("is_stable_closed_loop with zero gain on unstable system") {
    // Unstable system with zero gain -> still unstable
    Eigen::Matrix<double, 2, 2> A;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 2> K;
    A << 2.0, 0.0,
         0.0, 0.5;
    B << 1.0,
         0.0;
    K << 0.0, 0.0;

    CHECK_FALSE(ctrlpp::is_stable_closed_loop<Policy, double, 2, 1>(A, B, K));
}

TEST_CASE("is_stable_observer with stabilizing observer gain") {
    // A=[[0.9,0.1],[0,0.8]], C=[1,0], L chosen to stabilize
    // A-LC with L=[[0.5],[0.1]]:
    // A-LC = [[0.9-0.5, 0.1],[0-0.1, 0.8]] = [[0.4, 0.1],[-0.1, 0.8]]
    Eigen::Matrix<double, 2, 2> A;
    Eigen::Matrix<double, 1, 2> C;
    Eigen::Matrix<double, 2, 1> L;
    A << 0.9, 0.1,
         0.0, 0.8;
    C << 1.0, 0.0;
    L << 0.5,
         0.1;

    CHECK(ctrlpp::is_stable_observer<Policy, double, 2, 1>(A, L, C));
}

TEST_CASE("is_controllable with MIMO system") {
    // 3-state, 2-input system
    Eigen::Matrix<double, 3, 3> A;
    Eigen::Matrix<double, 3, 2> B;
    A << 0.9, 0.1, 0.0,
         0.0, 0.8, 0.2,
         0.0, 0.0, 0.7;
    B << 1.0, 0.0,
         0.0, 0.0,
         0.0, 1.0;

    // States 1 and 3 directly actuated; state 2 reachable through A coupling
    CHECK(ctrlpp::is_controllable<Policy, double, 3, 2>(A, B));
}
