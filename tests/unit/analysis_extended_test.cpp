#include "ctrlpp/model/analysis.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>

TEST_CASE("is_controllable with controllable double integrator")
{
    Eigen::Matrix<double, 2, 2> A;
    Eigen::Matrix<double, 2, 1> B;
    A << 1.0, 1.0, 0.0, 1.0;
    B << 0.5, 1.0;

    CHECK(ctrlpp::is_controllable<double, 2, 1>(A, B));
}

TEST_CASE("is_controllable with uncontrollable system")
{
    Eigen::Matrix<double, 2, 2> A;
    Eigen::Matrix<double, 2, 1> B;
    A << 2.0, 0.0, 0.0, 0.5;
    B << 0.0, 1.0;

    CHECK_FALSE(ctrlpp::is_controllable<double, 2, 1>(A, B));
}

TEST_CASE("is_observable with observable system")
{
    Eigen::Matrix<double, 2, 2> A;
    Eigen::Matrix<double, 1, 2> C;
    A << 1.0, 1.0, 0.0, 1.0;
    C << 1.0, 0.0;

    CHECK(ctrlpp::is_observable<double, 2, 1>(A, C));
}

TEST_CASE("is_observable with unobservable system")
{
    Eigen::Matrix<double, 2, 2> A;
    Eigen::Matrix<double, 1, 2> C;
    A << 2.0, 0.0, 0.0, 0.5;
    C << 0.0, 1.0;

    CHECK_FALSE(ctrlpp::is_observable<double, 2, 1>(A, C));
}

TEST_CASE("observability duality with controllability")
{
    Eigen::Matrix<double, 2, 2> A;
    Eigen::Matrix<double, 1, 2> C;
    A << 1.0, 1.0, 0.0, 1.0;
    C << 1.0, 0.0;

    auto At = A.transpose().eval();
    auto Ct = C.transpose().eval();

    bool obs = ctrlpp::is_observable<double, 2, 1>(A, C);
    bool ctrl_dual = ctrlpp::is_controllable<double, 2, 1>(At, Ct);
    CHECK(obs == ctrl_dual);
}

TEST_CASE("is_stable_closed_loop with stabilizing gain")
{
    Eigen::Matrix<double, 2, 2> A;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 2> K;
    A << 1.0, 1.0, 0.0, 1.0;
    B << 0.5, 1.0;
    K << 0.4, 0.8;

    Eigen::Matrix<double, 2, 2> Acl = A - B * K;
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> solver(Acl, false);
    bool all_stable = true;
    for(int i = 0; i < 2; ++i)
        if(std::abs(solver.eigenvalues()(i)) >= 1.0)
            all_stable = false;

    if(all_stable)
        CHECK(ctrlpp::is_stable_closed_loop<double, 2, 1>(A, B, K));
}

TEST_CASE("is_stable_closed_loop with zero gain on unstable system")
{
    Eigen::Matrix<double, 2, 2> A;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 2> K;
    A << 2.0, 0.0, 0.0, 0.5;
    B << 1.0, 0.0;
    K << 0.0, 0.0;

    CHECK_FALSE(ctrlpp::is_stable_closed_loop<double, 2, 1>(A, B, K));
}

TEST_CASE("is_stable_observer with stabilizing observer gain")
{
    Eigen::Matrix<double, 2, 2> A;
    Eigen::Matrix<double, 1, 2> C;
    Eigen::Matrix<double, 2, 1> L;
    A << 0.9, 0.1, 0.0, 0.8;
    C << 1.0, 0.0;
    L << 0.5, 0.1;

    CHECK(ctrlpp::is_stable_observer<double, 2, 1>(A, L, C));
}

TEST_CASE("is_controllable with MIMO system")
{
    Eigen::Matrix<double, 3, 3> A;
    Eigen::Matrix<double, 3, 2> B;
    A << 0.9, 0.1, 0.0, 0.0, 0.8, 0.2, 0.0, 0.0, 0.7;
    B << 1.0, 0.0, 0.0, 0.0, 0.0, 1.0;

    CHECK(ctrlpp::is_controllable<double, 3, 2>(A, B));
}
