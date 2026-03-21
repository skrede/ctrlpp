#include "ctrlpp/detail/numerical_diff.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace ctrlpp;
using namespace ctrlpp::detail;

namespace {

constexpr double tol = 1e-7;

// Linear dynamics: f(x, u) = A*x + B*u
struct LinearDynamics {
    Matrix<double, 2, 2> A;
    Matrix<double, 2, 1> B;

    auto operator()(const Vector<double, 2>& x,
                    const Vector<double, 1>& u) const -> Vector<double, 2> {
        return A * x + B * u;
    }
};

// Linear measurement: h(x) = C*x
struct LinearMeasurement {
    Matrix<double, 1, 2> C;

    auto operator()(const Vector<double, 2>& x) const -> Vector<double, 1> {
        return C * x;
    }
};

}

TEST_CASE("numerical_jacobian_x returns A for linear dynamics", "[numerical_diff]") {
    Matrix<double, 2, 2> A;
    A << 0.9, 0.1,
        -0.2, 0.8;
    Matrix<double, 2, 1> B;
    B << 1.0, 0.5;

    LinearDynamics dyn{A, B};
    Vector<double, 2> x;
    x << 1.0, -0.5;
    Vector<double, 1> u;
    u << 0.3;

    auto jac = numerical_jacobian_x<double, 2, 1>(dyn, x, u);

    for (Eigen::Index i = 0; i < 2; ++i) {
        for (Eigen::Index j = 0; j < 2; ++j) {
            REQUIRE_THAT(jac(i, j), Catch::Matchers::WithinAbs(A(i, j), tol));
        }
    }
}

TEST_CASE("numerical_jacobian_u returns B for linear dynamics", "[numerical_diff]") {
    Matrix<double, 2, 2> A;
    A << 0.9, 0.1,
        -0.2, 0.8;
    Matrix<double, 2, 1> B;
    B << 1.0, 0.5;

    LinearDynamics dyn{A, B};
    Vector<double, 2> x;
    x << 1.0, -0.5;
    Vector<double, 1> u;
    u << 0.3;

    auto jac = numerical_jacobian_u<double, 2, 1>(dyn, x, u);

    for (Eigen::Index i = 0; i < 2; ++i) {
        for (Eigen::Index j = 0; j < 1; ++j) {
            REQUIRE_THAT(jac(i, j), Catch::Matchers::WithinAbs(B(i, j), tol));
        }
    }
}

TEST_CASE("numerical_jacobian_h returns C for linear measurement", "[numerical_diff]") {
    Matrix<double, 1, 2> C;
    C << 0.7, -0.3;

    LinearMeasurement meas{C};
    Vector<double, 2> x;
    x << 2.0, 1.0;

    auto jac = numerical_jacobian_h<double, 2, 1>(meas, x);

    for (Eigen::Index i = 0; i < 1; ++i) {
        for (Eigen::Index j = 0; j < 2; ++j) {
            REQUIRE_THAT(jac(i, j), Catch::Matchers::WithinAbs(C(i, j), tol));
        }
    }
}

TEST_CASE("numerical_jacobian_x for nonlinear f(x) = x^2", "[numerical_diff]") {
    // f(x, u) = [x0^2, x1^2] -- Jacobian is diag(2*x0, 2*x1)
    auto nonlinear = [](const Vector<double, 2>& x,
                        const Vector<double, 1>&) -> Vector<double, 2> {
        Vector<double, 2> result;
        result[0] = x[0] * x[0];
        result[1] = x[1] * x[1];
        return result;
    };

    Vector<double, 2> x;
    x << 3.0, -2.0;
    Vector<double, 1> u;
    u << 0.0;

    auto jac = numerical_jacobian_x<double, 2, 1>(nonlinear, x, u);

    // Expected: diag(2*3, 2*(-2)) = diag(6, -4)
    REQUIRE_THAT(jac(0, 0), Catch::Matchers::WithinAbs(6.0, tol));
    REQUIRE_THAT(jac(0, 1), Catch::Matchers::WithinAbs(0.0, tol));
    REQUIRE_THAT(jac(1, 0), Catch::Matchers::WithinAbs(0.0, tol));
    REQUIRE_THAT(jac(1, 1), Catch::Matchers::WithinAbs(-4.0, tol));
}

TEST_CASE("numerical_jacobian_h for nonlinear h(x) = [sin(x0), cos(x1)]", "[numerical_diff]") {
    auto nonlinear_h = [](const Vector<double, 2>& x) -> Vector<double, 2> {
        Vector<double, 2> result;
        result[0] = std::sin(x[0]);
        result[1] = std::cos(x[1]);
        return result;
    };

    Vector<double, 2> x;
    x << 1.0, 0.5;

    auto jac = numerical_jacobian_h<double, 2, 2>(nonlinear_h, x);

    // Expected: [[cos(1), 0], [0, -sin(0.5)]]
    REQUIRE_THAT(jac(0, 0), Catch::Matchers::WithinAbs(std::cos(1.0), tol));
    REQUIRE_THAT(jac(0, 1), Catch::Matchers::WithinAbs(0.0, tol));
    REQUIRE_THAT(jac(1, 0), Catch::Matchers::WithinAbs(0.0, tol));
    REQUIRE_THAT(jac(1, 1), Catch::Matchers::WithinAbs(-std::sin(0.5), tol));
}
