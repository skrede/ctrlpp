// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./batch_arx_example' using 1:2 with lines title 'actual', '' using 1:3 with lines title 'predicted'"
// Redirect: ./batch_arx_example > output.csv

#include "ctrlpp/sysid.h"

#include <Eigen/Dense>

#include <random>
#include <cstddef>
#include <iostream>

int main()
{
    constexpr std::size_t N = 500;

    Eigen::Matrix<double, 1, static_cast<int>(N)> Y;
    Eigen::Matrix<double, 1, static_cast<int>(N)> U;

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> u_dist(-1.0, 1.0);

    double y = 0.0;
    double u_prev = 0.0;
    for(std::size_t t = 0; t < N; ++t)
    {
        double u = u_dist(gen);
        double y_new = 0.7 * y + 0.3 * u_prev;
        Y(0, static_cast<Eigen::Index>(t)) = y_new;
        U(0, static_cast<Eigen::Index>(t)) = u;
        y = y_new;
        u_prev = u;
    }

    // batch_arx refuses a record it cannot form a regressor from, rather than return a fit built on wrapped
    // index arithmetic.
    auto result = ctrlpp::batch_arx<1, 1>(Y, U);
    if(!result)
    {
        char const* reason = "unknown rejection";
        switch(result.error())
        {
            case ctrlpp::sysid_error::record_length_mismatch:
                reason = "the output and input records differ in length";
                break;
            case ctrlpp::sysid_error::record_not_single_row:
                reason = "a record is not a single row";
                break;
            case ctrlpp::sysid_error::too_few_samples:
                reason = "fewer samples than the model order requires";
                break;
            case ctrlpp::sysid_error::non_finite_sample:
                reason = "a sample is NaN or infinite";
                break;
            case ctrlpp::sysid_error::block_rows_below_order:
                reason = "too few block rows for the subspace model order";
                break;
            case ctrlpp::sysid_error::model_not_identified:
                reason = "the subspace identification produced no model";
                break;
        }
        std::cerr << "Batch ARX identification rejected the data: " << reason << '\n';
        return 1;
    }

    // A valid record can still fail to determine the model. That is not a refusal: the fit comes back with
    // the deficiency in its diagnostics, and the caller decides whether to use it.
    auto const& diagnostics = result->diagnostics;
    if(diagnostics.numerical_rank < diagnostics.parameter_count)
        std::cerr << "Batch ARX: the record determined only " << diagnostics.numerical_rank << " of " << diagnostics.parameter_count << " coefficients\n";

    auto const& sys = result->system;
    Eigen::Vector<double, 1> x = Eigen::Vector<double, 1>::Zero();

    std::cout << "# step,actual,predicted\n";
    for(std::size_t t = 1; t < N; ++t)
    {
        Eigen::Vector<double, 1> u_t;
        u_t(0) = U(0, static_cast<Eigen::Index>(t - 1));
        x = sys.A * x + sys.B * u_t;
        Eigen::Vector<double, 1> y_hat = sys.C * x + sys.D * u_t;
        std::cout << t << ',' << Y(0, static_cast<Eigen::Index>(t)) << ',' << y_hat(0) << '\n';
    }

    // stderr, so the model summary does not mix with the CSV.
    std::cerr << "Batch ARX: A=" << sys.A << " B=" << sys.B << " NRMSE=" << result->metrics.nrmse << " VAF=" << result->metrics.vaf << "%"
              << " fit residual=" << diagnostics.residual_norm << " over " << diagnostics.effective_samples << " rows\n";

    return 0;
}
