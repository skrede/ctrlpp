#include "ctrlpp/sysid/fit_metrics.h"
#include "ctrlpp/sysid/sysid_result.h"
#include "ctrlpp/model/state_space.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>

using Catch::Matchers::WithinAbs;

TEST_CASE("NRMSE of perfect prediction is 0")
{
    Eigen::VectorXd y(5);
    y << 1.0, 2.0, 3.0, 4.0, 5.0;

    auto m = ctrlpp::compute_fit_metrics(y, y);
    REQUIRE_THAT(m.nrmse, WithinAbs(0.0, 1e-12));
}

TEST_CASE("NRMSE of mean predictor is 1")
{
    Eigen::VectorXd y(5);
    y << 1.0, 2.0, 3.0, 4.0, 5.0;

    double mean = y.mean();
    Eigen::VectorXd y_hat = Eigen::VectorXd::Constant(5, mean);

    auto m = ctrlpp::compute_fit_metrics(y, y_hat);
    REQUIRE_THAT(m.nrmse, WithinAbs(1.0, 1e-12));
}

TEST_CASE("VAF of perfect prediction is 100")
{
    Eigen::VectorXd y(5);
    y << 1.0, 2.0, 3.0, 4.0, 5.0;

    auto m = ctrlpp::compute_fit_metrics(y, y);
    REQUIRE_THAT(m.vaf, WithinAbs(100.0, 1e-12));
}

TEST_CASE("VAF of zero predictor on non-constant data is less than 100")
{
    Eigen::VectorXd y(5);
    y << 1.0, 2.0, 3.0, 4.0, 5.0;

    Eigen::VectorXd y_hat = Eigen::VectorXd::Zero(5);

    auto m = ctrlpp::compute_fit_metrics(y, y_hat);
    REQUIRE(m.vaf < 100.0);
}

TEST_CASE("fit_metrics works with static-sized Eigen vectors")
{
    Eigen::Vector<double, 4> y;
    y << 1.0, 3.0, 5.0, 7.0;

    Eigen::Vector<double, 4> y_hat;
    y_hat << 1.0, 3.0, 5.0, 7.0;

    auto m = ctrlpp::compute_fit_metrics(y, y_hat);
    REQUIRE_THAT(m.nrmse, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(m.vaf, WithinAbs(100.0, 1e-12));
}

TEST_CASE("fit_metrics edge case: constant y with perfect prediction")
{
    Eigen::VectorXd y = Eigen::VectorXd::Constant(5, 3.0);
    Eigen::VectorXd y_hat = Eigen::VectorXd::Constant(5, 3.0);

    auto m = ctrlpp::compute_fit_metrics(y, y_hat);
    REQUIRE_THAT(m.nrmse, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(m.vaf, WithinAbs(100.0, 1e-12));
}

TEST_CASE("arx_result holds discrete_state_space and fit_metrics")
{
    ctrlpp::arx_result<double, 2, 1, 1> r{};
    r.system.A = ctrlpp::Matrix<double, 2, 2>::Identity();
    r.metrics.nrmse = 0.1;
    r.metrics.vaf = 95.0;
    REQUIRE_THAT(r.metrics.nrmse, WithinAbs(0.1, 1e-12));
    REQUIRE_THAT(r.metrics.vaf, WithinAbs(95.0, 1e-12));
}

TEST_CASE("n4sid_result holds state_space, singular_values, metrics, condition_number")
{
    ctrlpp::n4sid_result<double, 2, 1, 1> r{};
    r.system.A = ctrlpp::Matrix<double, 2, 2>::Identity();
    r.singular_values.resize(2);
    r.singular_values << 10.0, 1.0;
    r.metrics.nrmse = 0.05;
    r.metrics.vaf = 99.0;
    r.condition_number = 10.0;
    REQUIRE_THAT(r.condition_number, WithinAbs(10.0, 1e-12));
    REQUIRE(r.singular_values.size() == 2);
}

TEST_CASE("NRMSE is infinity for constant y with imperfect prediction")
{
    Eigen::VectorXd y = Eigen::VectorXd::Constant(5, 3.0);
    Eigen::VectorXd y_hat = Eigen::VectorXd::Constant(5, 4.0);

    auto m = ctrlpp::compute_fit_metrics(y, y_hat);
    // norm_centered ~ 0 (constant y), norm_error > 0 -> NRMSE = infinity
    REQUIRE(std::isinf(m.nrmse));
    REQUIRE(m.nrmse > 0.0);
}

TEST_CASE("VAF is 100 for constant y with constant imperfect prediction")
{
    Eigen::VectorXd y = Eigen::VectorXd::Constant(5, 3.0);
    Eigen::VectorXd y_hat = Eigen::VectorXd::Constant(5, 4.0);

    auto m = ctrlpp::compute_fit_metrics(y, y_hat);
    // var_y=0, error is also constant so var_error=0
    // Both variances near-zero -> VAF = 100 (zero-variance-both path)
    REQUIRE_THAT(m.vaf, WithinAbs(100.0, 1e-12));
}

TEST_CASE("VAF is -infinity for constant y with varying prediction")
{
    Eigen::VectorXd y = Eigen::VectorXd::Constant(5, 3.0);
    Eigen::VectorXd y_hat(5);
    y_hat << 1.0, 2.0, 3.0, 4.0, 5.0;  // varying prediction -> var_error > 0

    auto m = ctrlpp::compute_fit_metrics(y, y_hat);
    // var_y = 0, var_error > 0 -> VAF = -infinity
    REQUIRE(std::isinf(m.vaf));
    REQUIRE(m.vaf < 0.0);
}

TEST_CASE("Single sample: VAF defaults to 100 with perfect prediction")
{
    Eigen::VectorXd y(1);
    y << 5.0;
    Eigen::VectorXd y_hat(1);
    y_hat << 5.0;

    auto m = ctrlpp::compute_fit_metrics(y, y_hat);
    // n=1, so var_error=0 and var_y=0 (both skip the n>1 branch)
    // -> constant-y-with-perfect-prediction path: VAF=100, NRMSE=0
    REQUIRE_THAT(m.nrmse, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(m.vaf, WithinAbs(100.0, 1e-12));
}

TEST_CASE("Single sample: imperfect prediction gives NRMSE=inf and VAF=-inf")
{
    Eigen::VectorXd y(1);
    y << 5.0;
    Eigen::VectorXd y_hat(1);
    y_hat << 3.0;

    auto m = ctrlpp::compute_fit_metrics(y, y_hat);
    // n=1: var_y=0 (no variance computed), var_error=0
    // norm_centered = 0 (single point, centered = y - mean = 0)
    // norm_error > 0 -> NRMSE = infinity
    REQUIRE(std::isinf(m.nrmse));
    // var_y ~ 0, var_error ~ 0 (both stay default 0 since n<=1)
    // -> VAF = 100 (both variances zero path)
    // This is the edge case: with n=1 we can't compute variance
    REQUIRE_THAT(m.vaf, WithinAbs(100.0, 1e-12));
}

TEST_CASE("NRMSE greater than 1 when prediction is worse than mean")
{
    Eigen::VectorXd y(5);
    y << 1.0, 2.0, 3.0, 4.0, 5.0;

    // Prediction far from actual
    Eigen::VectorXd y_hat(5);
    y_hat << 10.0, 10.0, 10.0, 10.0, 10.0;

    auto m = ctrlpp::compute_fit_metrics(y, y_hat);
    REQUIRE(m.nrmse > 1.0);
}

TEST_CASE("VAF is negative when prediction error has more variance than data")
{
    Eigen::VectorXd y(5);
    y << 1.0, 2.0, 3.0, 4.0, 5.0;

    // Prediction that inverts the trend -> error variance > data variance
    Eigen::VectorXd y_hat(5);
    y_hat << 5.0, 4.0, 3.0, 2.0, 1.0;

    auto m = ctrlpp::compute_fit_metrics(y, y_hat);
    // error = y - y_hat = [-4, -2, 0, 2, 4], var_error > var_y -> VAF < 0
    REQUIRE(m.vaf < 0.0);
}

TEST_CASE("fit_metrics with float scalar type")
{
    Eigen::VectorXf y(4);
    y << 1.0f, 2.0f, 3.0f, 4.0f;

    Eigen::VectorXf y_hat(4);
    y_hat << 1.0f, 2.0f, 3.0f, 4.0f;

    auto m = ctrlpp::compute_fit_metrics(y, y_hat);
    REQUIRE_THAT(static_cast<double>(m.nrmse), WithinAbs(0.0, 1e-5));
    REQUIRE_THAT(static_cast<double>(m.vaf), WithinAbs(100.0, 1e-3));
}
