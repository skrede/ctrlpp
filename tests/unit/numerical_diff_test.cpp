#include "ctrlpp/detail/numerical_diff.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using namespace ctrlpp;
using namespace ctrlpp::detail;

namespace
{

constexpr double tol = 1e-7;

// Linear dynamics: f(x, u) = A*x + B*u
struct LinearDynamics
{
    Matrix<double, 2, 2> A;
    Matrix<double, 2, 1> B;

    auto operator()(const Vector<double, 2>& x, const Vector<double, 1>& u) const -> Vector<double, 2> { return A * x + B * u; }
};

// Linear measurement: h(x) = C*x
struct LinearMeasurement
{
    Matrix<double, 1, 2> C;

    auto operator()(const Vector<double, 2>& x) const -> Vector<double, 1> { return C * x; }
};

} // namespace

TEST_CASE("numerical_jacobian_x returns A for linear dynamics", "[numerical_diff]")
{
    Matrix<double, 2, 2> A;
    A << 0.9, 0.1, -0.2, 0.8;
    Matrix<double, 2, 1> B;
    B << 1.0, 0.5;

    LinearDynamics dyn{A, B};
    Vector<double, 2> x;
    x << 1.0, -0.5;
    Vector<double, 1> u;
    u << 0.3;

    auto jac = numerical_jacobian_x<double, 2, 1>(dyn, x, u);

    for(Eigen::Index i = 0; i < 2; ++i)
    {
        for(Eigen::Index j = 0; j < 2; ++j)
        {
            REQUIRE_THAT(jac(i, j), Catch::Matchers::WithinAbs(A(i, j), tol));
        }
    }
}

TEST_CASE("numerical_jacobian_u returns B for linear dynamics", "[numerical_diff]")
{
    Matrix<double, 2, 2> A;
    A << 0.9, 0.1, -0.2, 0.8;
    Matrix<double, 2, 1> B;
    B << 1.0, 0.5;

    LinearDynamics dyn{A, B};
    Vector<double, 2> x;
    x << 1.0, -0.5;
    Vector<double, 1> u;
    u << 0.3;

    auto jac = numerical_jacobian_u<double, 2, 1>(dyn, x, u);

    for(Eigen::Index i = 0; i < 2; ++i)
    {
        for(Eigen::Index j = 0; j < 1; ++j)
        {
            REQUIRE_THAT(jac(i, j), Catch::Matchers::WithinAbs(B(i, j), tol));
        }
    }
}

TEST_CASE("numerical_jacobian_h returns C for linear measurement", "[numerical_diff]")
{
    Matrix<double, 1, 2> C;
    C << 0.7, -0.3;

    LinearMeasurement meas{C};
    Vector<double, 2> x;
    x << 2.0, 1.0;

    auto jac = numerical_jacobian_h<double, 2, 1>(meas, x);

    for(Eigen::Index i = 0; i < 1; ++i)
    {
        for(Eigen::Index j = 0; j < 2; ++j)
        {
            REQUIRE_THAT(jac(i, j), Catch::Matchers::WithinAbs(C(i, j), tol));
        }
    }
}

TEST_CASE("numerical_jacobian_x for nonlinear f(x) = x^2", "[numerical_diff]")
{
    // f(x, u) = [x0^2, x1^2] -- Jacobian is diag(2*x0, 2*x1)
    auto nonlinear = [](const Vector<double, 2>& x, const Vector<double, 1>&) -> Vector<double, 2>
    {
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

TEST_CASE("numerical_jacobian_h for nonlinear h(x) = [sin(x0), cos(x1)]", "[numerical_diff]")
{
    auto nonlinear_h = [](const Vector<double, 2>& x) -> Vector<double, 2>
    {
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

// ---------------------------------------------------------------------------
// Tests for previously untested functions
// ---------------------------------------------------------------------------

TEST_CASE("finite_diff_gradient for quadratic scalar function", "[numerical_diff]")
{
    // f(z) = z0^2 + 3*z1 => grad = [2*z0, 3]
    std::function<double(std::span<const double>)> f = [](std::span<const double> z) -> double
    {
        return z[0] * z[0] + 3.0 * z[1];
    };

    std::vector<double> z = {2.0, 5.0};
    std::vector<double> grad(2, 0.0);
    finite_diff_gradient<double>(f, z, grad);

    REQUIRE_THAT(grad[0], Catch::Matchers::WithinAbs(4.0, tol));
    REQUIRE_THAT(grad[1], Catch::Matchers::WithinAbs(3.0, tol));
}

TEST_CASE("finite_diff_gradient at zero", "[numerical_diff]")
{
    // f(z) = z0^2 + z1^2 => grad at origin = [0, 0]
    std::function<double(std::span<const double>)> f = [](std::span<const double> z) -> double
    {
        return z[0] * z[0] + z[1] * z[1];
    };

    std::vector<double> z = {0.0, 0.0};
    std::vector<double> grad(2, 0.0);
    finite_diff_gradient<double>(f, z, grad);

    REQUIRE_THAT(grad[0], Catch::Matchers::WithinAbs(0.0, tol));
    REQUIRE_THAT(grad[1], Catch::Matchers::WithinAbs(0.0, tol));
}

TEST_CASE("finite_diff_gradient with large values", "[numerical_diff]")
{
    // f(z) = z0 * z1 => grad = [z1, z0]
    std::function<double(std::span<const double>)> f = [](std::span<const double> z) -> double
    {
        return z[0] * z[1];
    };

    std::vector<double> z = {1e6, 1e6};
    std::vector<double> grad(2, 0.0);
    finite_diff_gradient<double>(f, z, grad);

    REQUIRE_THAT(grad[0], Catch::Matchers::WithinAbs(1e6, 1e-1));
    REQUIRE_THAT(grad[1], Catch::Matchers::WithinAbs(1e6, 1e-1));
}

TEST_CASE("finite_diff_jacobian for linear constraint c(z) = [2*z0 + z1, z0 - z1]", "[numerical_diff]")
{
    // Jacobian should be [[2, 1], [1, -1]]
    std::function<void(std::span<const double>, std::span<double>)> c =
        [](std::span<const double> z, std::span<double> out)
    {
        out[0] = 2.0 * z[0] + z[1];
        out[1] = z[0] - z[1];
    };

    std::vector<double> z = {1.0, 2.0};
    std::vector<double> jac(4, 0.0);
    finite_diff_jacobian<double>(c, 2, z, jac);

    // Row-major: jac[i*n + j]
    REQUIRE_THAT(jac[0], Catch::Matchers::WithinAbs(2.0, tol)); // dc0/dz0
    REQUIRE_THAT(jac[1], Catch::Matchers::WithinAbs(1.0, tol)); // dc0/dz1
    REQUIRE_THAT(jac[2], Catch::Matchers::WithinAbs(1.0, tol)); // dc1/dz0
    REQUIRE_THAT(jac[3], Catch::Matchers::WithinAbs(-1.0, tol)); // dc1/dz1
}

TEST_CASE("finite_diff_jacobian at zero point", "[numerical_diff]")
{
    // c(z) = [z0*z1, z0^2] => Jac at (0,0) = [[0, 0], [0, 0]]
    std::function<void(std::span<const double>, std::span<double>)> c =
        [](std::span<const double> z, std::span<double> out)
    {
        out[0] = z[0] * z[1];
        out[1] = z[0] * z[0];
    };

    std::vector<double> z = {0.0, 0.0};
    std::vector<double> jac(4, 0.0);
    finite_diff_jacobian<double>(c, 2, z, jac);

    for(int i = 0; i < 4; ++i)
    {
        REQUIRE_THAT(jac[i], Catch::Matchers::WithinAbs(0.0, tol));
    }
}

TEST_CASE("numerical_jacobian_g for linear constraint g(x, u) = [x0 + u0, x1 - u0]", "[numerical_diff]")
{
    // NX=2, NU=1, NC=2. Jacobian = [[1, 0, 1], [0, 1, -1]]
    auto g = [](const Vector<double, 2>& x, const Vector<double, 1>& u) -> Vector<double, 2>
    {
        Vector<double, 2> result;
        result[0] = x[0] + u[0];
        result[1] = x[1] - u[0];
        return result;
    };

    Vector<double, 2> x;
    x << 1.0, 2.0;
    Vector<double, 1> u;
    u << 0.5;

    auto jac = numerical_jacobian_g<double, 2, 1, 2>(g, x, u);

    // 2 x 3 matrix
    REQUIRE_THAT(jac(0, 0), Catch::Matchers::WithinAbs(1.0, tol)); // dg0/dx0
    REQUIRE_THAT(jac(0, 1), Catch::Matchers::WithinAbs(0.0, tol)); // dg0/dx1
    REQUIRE_THAT(jac(0, 2), Catch::Matchers::WithinAbs(1.0, tol)); // dg0/du0
    REQUIRE_THAT(jac(1, 0), Catch::Matchers::WithinAbs(0.0, tol)); // dg1/dx0
    REQUIRE_THAT(jac(1, 1), Catch::Matchers::WithinAbs(1.0, tol)); // dg1/dx1
    REQUIRE_THAT(jac(1, 2), Catch::Matchers::WithinAbs(-1.0, tol)); // dg1/du0
}

TEST_CASE("numerical_jacobian_g for nonlinear constraint", "[numerical_diff]")
{
    // g(x, u) = [x0*x1 + u0^2], NX=2, NU=1, NC=1
    // dg/dx0 = x1, dg/dx1 = x0, dg/du0 = 2*u0
    auto g = [](const Vector<double, 2>& x, const Vector<double, 1>& u) -> Vector<double, 1>
    {
        Vector<double, 1> result;
        result[0] = x[0] * x[1] + u[0] * u[0];
        return result;
    };

    Vector<double, 2> x;
    x << 3.0, -2.0;
    Vector<double, 1> u;
    u << 1.5;

    auto jac = numerical_jacobian_g<double, 2, 1, 1>(g, x, u);

    REQUIRE_THAT(jac(0, 0), Catch::Matchers::WithinAbs(-2.0, tol)); // dg/dx0 = x1
    REQUIRE_THAT(jac(0, 1), Catch::Matchers::WithinAbs(3.0, tol)); // dg/dx1 = x0
    REQUIRE_THAT(jac(0, 2), Catch::Matchers::WithinAbs(3.0, tol)); // dg/du0 = 2*u0
}

TEST_CASE("numerical_jacobian_tc delegates to numerical_jacobian_h", "[numerical_diff]")
{
    // Terminal constraint: h(x) = [x0 + x1], NX=2, NTC=1
    auto h = [](const Vector<double, 2>& x) -> Vector<double, 1>
    {
        Vector<double, 1> result;
        result[0] = x[0] + x[1];
        return result;
    };

    Vector<double, 2> x;
    x << 1.0, 2.0;

    auto jac_tc = numerical_jacobian_tc<double, 2, 1>(h, x);
    auto jac_h = numerical_jacobian_h<double, 2, 1>(h, x);

    for(Eigen::Index i = 0; i < 1; ++i)
    {
        for(Eigen::Index j = 0; j < 2; ++j)
        {
            REQUIRE_THAT(jac_tc(i, j), Catch::Matchers::WithinAbs(jac_h(i, j), 1e-15));
        }
    }
}

TEST_CASE("numerical_jacobian_x at zero state", "[numerical_diff]")
{
    // f(x, u) = [x0 + x1, x0*x1] => dF/dx at (0,0) = [[1, 1], [0, 0]]
    auto f = [](const Vector<double, 2>& x, const Vector<double, 1>&) -> Vector<double, 2>
    {
        Vector<double, 2> result;
        result[0] = x[0] + x[1];
        result[1] = x[0] * x[1];
        return result;
    };

    Vector<double, 2> x = Vector<double, 2>::Zero();
    Vector<double, 1> u = Vector<double, 1>::Zero();

    auto jac = numerical_jacobian_x<double, 2, 1>(f, x, u);

    REQUIRE_THAT(jac(0, 0), Catch::Matchers::WithinAbs(1.0, tol));
    REQUIRE_THAT(jac(0, 1), Catch::Matchers::WithinAbs(1.0, tol));
    REQUIRE_THAT(jac(1, 0), Catch::Matchers::WithinAbs(0.0, tol));
    REQUIRE_THAT(jac(1, 1), Catch::Matchers::WithinAbs(0.0, tol));
}

TEST_CASE("numerical_jacobian_h with multi-output measurement", "[numerical_diff]")
{
    // h(x) = [x0, x1, x0+x1] for NX=2, NY=3
    auto h = [](const Vector<double, 2>& x) -> Vector<double, 3>
    {
        Vector<double, 3> result;
        result[0] = x[0];
        result[1] = x[1];
        result[2] = x[0] + x[1];
        return result;
    };

    Vector<double, 2> x;
    x << 5.0, -3.0;

    auto jac = numerical_jacobian_h<double, 2, 3>(h, x);

    // Expected: [[1, 0], [0, 1], [1, 1]]
    REQUIRE_THAT(jac(0, 0), Catch::Matchers::WithinAbs(1.0, tol));
    REQUIRE_THAT(jac(0, 1), Catch::Matchers::WithinAbs(0.0, tol));
    REQUIRE_THAT(jac(1, 0), Catch::Matchers::WithinAbs(0.0, tol));
    REQUIRE_THAT(jac(1, 1), Catch::Matchers::WithinAbs(1.0, tol));
    REQUIRE_THAT(jac(2, 0), Catch::Matchers::WithinAbs(1.0, tol));
    REQUIRE_THAT(jac(2, 1), Catch::Matchers::WithinAbs(1.0, tol));
}

TEST_CASE("numerical_jacobian_u with multi-input", "[numerical_diff]")
{
    // f(x, u) = [x0 + u0 - u1, x1 + 2*u1], NX=2, NU=2
    // dF/du = [[1, -1], [0, 2]]
    auto f = [](const Vector<double, 2>& x, const Vector<double, 2>& u) -> Vector<double, 2>
    {
        Vector<double, 2> result;
        result[0] = x[0] + u[0] - u[1];
        result[1] = x[1] + 2.0 * u[1];
        return result;
    };

    Vector<double, 2> x;
    x << 1.0, 2.0;
    Vector<double, 2> u;
    u << 0.5, -0.3;

    auto jac = numerical_jacobian_u<double, 2, 2>(f, x, u);

    REQUIRE_THAT(jac(0, 0), Catch::Matchers::WithinAbs(1.0, tol));
    REQUIRE_THAT(jac(0, 1), Catch::Matchers::WithinAbs(-1.0, tol));
    REQUIRE_THAT(jac(1, 0), Catch::Matchers::WithinAbs(0.0, tol));
    REQUIRE_THAT(jac(1, 1), Catch::Matchers::WithinAbs(2.0, tol));
}
