#include "ctrlpp/traj/detail/tridiagonal.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>
#include <cmath>
#include <vector>

using Catch::Matchers::WithinAbs;

// -- Tridiagonal solver tests -------------------------------------------------

TEST_CASE("thomas_solve: 3x3 system", "[traj][tridiagonal]")
{
    // System: [[2,1,0],[1,3,1],[0,1,2]] * x = [1,2,3]
    // Dense solve: x = [0.1, 0.8, 1.1]  (verify by hand)
    std::vector<double> a = {0.0, 1.0, 1.0};    // sub-diagonal (a[0] unused)
    std::vector<double> b = {2.0, 3.0, 2.0};    // main diagonal
    std::vector<double> c = {1.0, 1.0, 0.0};    // super-diagonal (c[n-1] unused)
    std::vector<double> d = {1.0, 2.0, 3.0};    // RHS

    // Verify against Eigen dense solve
    Eigen::Matrix3d A;
    A << 2, 1, 0,
         1, 3, 1,
         0, 1, 2;
    Eigen::Vector3d rhs(1.0, 2.0, 3.0);
    Eigen::Vector3d x_ref = A.colPivHouseholderQr().solve(rhs);

    ctrlpp::detail::thomas_solve(a, b, c, d);

    for (int i = 0; i < 3; ++i) {
        REQUIRE_THAT(d[static_cast<std::size_t>(i)], WithinAbs(x_ref(i), 1e-12));
    }
}

TEST_CASE("thomas_solve: 5x5 diagonal-dominant system", "[traj][tridiagonal]")
{
    // Build a 5x5 diagonally dominant tridiagonal system
    constexpr int n = 5;
    std::vector<double> a(n), b(n), c(n), d(n);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd rhs(n);

    for (int i = 0; i < n; ++i) {
        b[static_cast<std::size_t>(i)] = 4.0;
        A(i, i) = 4.0;
        if (i > 0) {
            a[static_cast<std::size_t>(i)] = 1.0;
            A(i, i - 1) = 1.0;
        }
        if (i < n - 1) {
            c[static_cast<std::size_t>(i)] = 1.0;
            A(i, i + 1) = 1.0;
        }
        d[static_cast<std::size_t>(i)] = static_cast<double>(i + 1);
        rhs(i) = static_cast<double>(i + 1);
    }

    Eigen::VectorXd x_ref = A.colPivHouseholderQr().solve(rhs);

    ctrlpp::detail::thomas_solve(a, b, c, d);

    for (int i = 0; i < n; ++i) {
        REQUIRE_THAT(d[static_cast<std::size_t>(i)], WithinAbs(x_ref(i), 1e-12));
    }
}

TEST_CASE("cyclic_thomas_solve: 3x3 cyclic system", "[traj][tridiagonal]")
{
    // 3x3 cyclic tridiagonal:
    // [[4, 1, 0.5],
    //  [1, 4, 1  ],
    //  [0.7, 1, 4]]
    // Corner: A[0][n-1] = alpha = 0.5, A[n-1][0] = beta = 0.7
    constexpr int n = 3;
    std::vector<double> a = {0.0, 1.0, 1.0};
    std::vector<double> b = {4.0, 4.0, 4.0};
    std::vector<double> c = {1.0, 1.0, 0.0};
    std::vector<double> d = {1.0, 2.0, 3.0};
    double alpha = 0.5;  // A[0][n-1]
    double beta = 0.7;   // A[n-1][0]

    Eigen::Matrix3d A;
    A << 4.0, 1.0, 0.5,
         1.0, 4.0, 1.0,
         0.7, 1.0, 4.0;
    Eigen::Vector3d rhs(1.0, 2.0, 3.0);
    Eigen::Vector3d x_ref = A.colPivHouseholderQr().solve(rhs);

    ctrlpp::detail::cyclic_thomas_solve(a, b, c, d, alpha, beta);

    for (int i = 0; i < n; ++i) {
        REQUIRE_THAT(d[static_cast<std::size_t>(i)], WithinAbs(x_ref(i), 1e-12));
    }
}

TEST_CASE("cyclic_thomas_solve: 5x5 cyclic system", "[traj][tridiagonal]")
{
    constexpr int n = 5;
    std::vector<double> a(n), b(n), c(n), d(n);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd rhs(n);

    for (int i = 0; i < n; ++i) {
        b[static_cast<std::size_t>(i)] = 6.0;
        A(i, i) = 6.0;
        if (i > 0) {
            a[static_cast<std::size_t>(i)] = 1.0;
            A(i, i - 1) = 1.0;
        }
        if (i < n - 1) {
            c[static_cast<std::size_t>(i)] = 1.0;
            A(i, i + 1) = 1.0;
        }
        d[static_cast<std::size_t>(i)] = static_cast<double>(i * i + 1);
        rhs(i) = static_cast<double>(i * i + 1);
    }

    double alpha = 0.8;  // A[0][n-1]
    double beta = 0.3;   // A[n-1][0]
    A(0, n - 1) = alpha;
    A(n - 1, 0) = beta;

    Eigen::VectorXd x_ref = A.colPivHouseholderQr().solve(rhs);

    ctrlpp::detail::cyclic_thomas_solve(a, b, c, d, alpha, beta);

    for (int i = 0; i < n; ++i) {
        REQUIRE_THAT(d[static_cast<std::size_t>(i)], WithinAbs(x_ref(i), 1e-12));
    }
}
