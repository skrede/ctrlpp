#include "ctrlpp/control/care.h"
#include "ctrlpp/detail/care_methods.h"
#include "ctrlpp/detail/schur_reorder.h"
#include "ctrlpp/detail/hamiltonian_balance.h"


#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <tuple>
#include <cmath>
#include <limits>
#include <complex>


using care_method_tags = std::tuple<
    ctrlpp::detail::schur_care_method,
    ctrlpp::detail::sign_function_care_method,
    ctrlpp::detail::balanced_schur_care_method>;


TEMPLATE_LIST_TEST_CASE("care scalar integrator analytical", "[care]", care_method_tags)
{
    // A=0, B=1, Q=1, R=1 => P^2 = 1 => P = 1
    Eigen::Matrix<double, 1, 1> A, B, Q, R;
    A(0, 0) = 0.0;
    B(0, 0) = 1.0;
    Q(0, 0) = 1.0;
    R(0, 0) = 1.0;

    auto result = ctrlpp::care<double, 1, 1, TestType>(A, B, Q, R);
    REQUIRE(result.has_value());

    CHECK_THAT(result->P(0, 0), Catch::Matchers::WithinAbs(1.0, 1e-10));
}

TEMPLATE_LIST_TEST_CASE("care double integrator analytical", "[care]", care_method_tags)
{
    // Continuous double integrator: A=[[0,1],[0,0]], B=[[0],[1]], Q=I, R=1
    // Analytical P = [[sqrt(3), 1], [1, sqrt(3)]]
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 0.0, 1.0, 0.0, 0.0;
    B << 0.0, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;

    auto result = ctrlpp::care<double, 2, 1, TestType>(A, B, Q, R);
    REQUIRE(result.has_value());

    const auto& P = result->P;
    const double s3 = std::sqrt(3.0);

    CHECK_THAT(P(0, 0), Catch::Matchers::WithinAbs(s3, 1e-10));
    CHECK_THAT(P(1, 1), Catch::Matchers::WithinAbs(s3, 1e-10));
    CHECK_THAT(P(0, 1), Catch::Matchers::WithinAbs(1.0, 1e-10));
    CHECK_THAT(P(1, 0), Catch::Matchers::WithinAbs(1.0, 1e-10));

    auto residual = A.transpose() * P + P * A - P * B * R.inverse() * B.transpose() * P + Q;
    CHECK(residual.norm() < 1e-10);
}

TEMPLATE_LIST_TEST_CASE("care 3-state damped system", "[care]", care_method_tags)
{
    Eigen::Matrix<double, 3, 3> A, Q;
    Eigen::Matrix<double, 3, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << -0.5, 1.0, 0.0,
         0.0, -0.8, 1.0,
         0.0, 0.0, -1.0;
    B << 0.0, 0.0, 1.0;
    Q = Eigen::Matrix<double, 3, 3>::Identity();
    R(0, 0) = 1.0;

    auto result = ctrlpp::care<double, 3, 1, TestType>(A, B, Q, R);
    REQUIRE(result.has_value());

    const auto& P = result->P;

    CHECK((P - P.transpose()).norm() < 1e-10);

    auto residual = A.transpose() * P + P * A - P * B * R.inverse() * B.transpose() * P + Q;
    CHECK(residual.norm() < 1e-10);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> eigsolver(P);
    for(int i = 0; i < 3; ++i)
        CHECK(eigsolver.eigenvalues()(i) >= -1e-10);
}

TEMPLATE_LIST_TEST_CASE("care closed-loop stable", "[care]", care_method_tags)
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 0.0, 1.0, -0.1, -0.2;
    B << 0.0, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;

    auto result = ctrlpp::care<double, 2, 1, TestType>(A, B, Q, R);
    REQUIRE(result.has_value());

    const auto& P = result->P;
    Eigen::Matrix<double, 1, 2> K = R.inverse() * B.transpose() * P;
    Eigen::Matrix<double, 2, 2> Acl = A - B * K;

    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> solver(Acl, false);
    for(int i = 0; i < 2; ++i)
        CHECK(solver.eigenvalues()(i).real() < 0.0);
}

TEMPLATE_LIST_TEST_CASE("care non-stabilizable returns nullopt", "[care]", care_method_tags)
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 2.0, 0.0, 0.0, -0.5;
    B << 0.0, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;

    auto result = ctrlpp::care<double, 2, 1, TestType>(A, B, Q, R);
    CHECK_FALSE(result.has_value());
}

TEMPLATE_LIST_TEST_CASE("care with N cross-weight", "[care]", care_method_tags)
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B, N;
    Eigen::Matrix<double, 1, 1> R;

    A << 0.0, 1.0, -0.5, -0.3;
    B << 0.0, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;
    N << 0.1, 0.2;

    auto result_with_n = ctrlpp::care<double, 2, 1, TestType>(A, B, Q, R, N);
    REQUIRE(result_with_n.has_value());

    auto Rinv = R.inverse();
    Eigen::Matrix<double, 2, 2> Qp = Q - N * Rinv * N.transpose();
    Eigen::Matrix<double, 2, 2> Ap = A - B * Rinv * N.transpose();

    auto result_standard = ctrlpp::care<double, 2, 1, TestType>(Ap, B, Qp, R);
    REQUIRE(result_standard.has_value());

    CHECK((result_with_n->P - result_standard->P).norm() < 1e-10);
}

TEMPLATE_LIST_TEST_CASE("care solution is symmetric", "[care]", care_method_tags)
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << -0.5, 1.0, 0.0, -0.8;
    B << 1.0, 0.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R(0, 0) = 1.0;

    auto result = ctrlpp::care<double, 2, 1, TestType>(A, B, Q, R);
    REQUIRE(result.has_value());

    CHECK((result->P - result->P.transpose()).norm() < 1e-10);
}

TEST_CASE("care balanced_schur preserves subspace residual")
{
    using Scalar = double;
    constexpr int n  = 3;
    constexpr int n2 = 2 * n;
    using Mat_n      = Eigen::Matrix<Scalar, n, n>;
    using Mat_nu     = Eigen::Matrix<Scalar, n, 1>;
    using Mat_r      = Eigen::Matrix<Scalar, 1, 1>;
    using Mat_2n     = Eigen::Matrix<Scalar, n2, n2>;
    using Vec_2n     = Eigen::Matrix<Scalar, n2, 1>;

    Mat_n A;
    A << -0.5, 1.0, 0.0,
          0.0, -0.8, 1.0,
          0.0, 0.0, -1.0;
    Mat_nu B; B << 0.0, 0.0, 1.0;
    Mat_n  Q = Mat_n::Identity();
    Mat_r  R; R(0, 0) = 1.0;

    auto H_result = ctrlpp::detail::build_care_hamiltonian<Scalar, 3, 1>(A, B, Q, R);
    REQUIRE(H_result.has_value());
    const Mat_2n H_original = *H_result;

    Mat_2n H = H_original;
    Vec_2n D;
    ctrlpp::detail::balance_hamiltonian<Scalar, 3>(H, D);

    Eigen::RealSchur<Mat_2n> schur(H);
    REQUIRE(schur.info() == Eigen::Success);
    Mat_2n T = schur.matrixT();
    Mat_2n U = schur.matrixU();

    const Scalar scale      = T.cwiseAbs().maxCoeff();
    const Scalar eps        = std::numeric_limits<Scalar>::epsilon();
    const Scalar lhp_margin = eps * std::max(Scalar{1}, scale);
    auto predicate = [lhp_margin](std::complex<Scalar> lam) -> bool
    {
        return lam.real() < -lhp_margin;
    };
    auto rr = ctrlpp::detail::reorder_real_schur<Scalar, n2>(T, U, predicate);
    REQUIRE(rr.placed >= n);

    U.leftCols(n).array().colwise() *= D.array();

    const Scalar residual = (H_original * U.leftCols(n)
                             - U.leftCols(n) * T.topLeftCorner(n, n)).norm();
    const Scalar bound    = eps * Scalar{n2} * H_original.norm();
    CHECK(residual <= bound);
}
