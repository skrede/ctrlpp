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
