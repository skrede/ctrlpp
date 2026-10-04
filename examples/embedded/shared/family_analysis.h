#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_FAMILY_ANALYSIS_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_FAMILY_ANALYSIS_H

#include "dsp_demo.h"
#include "family_bounds.h"
#include "estimation_demo.h"
#include "trajectory_demo.h"
#include "golden_reference.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <algorithm>

namespace ctrlpp {

struct family_analysis
{
    double reference;
    double scale;
    const char *measured_name;
    double measured;
    const char *premise_name;
    double premise;
    double premise_limit;
};

namespace detail {

using analysis_matrix = Eigen::Matrix2d;
using analysis_maps   = std::array<analysis_matrix, kGoldenSteps>;

inline double infinity_norm(const analysis_matrix &m)
{
    return m.cwiseAbs().rowwise().sum().maxCoeff();
}

// The first-order map of a covariance error is X -> M*X*M^T, whose gain in the
// largest-entry norm is at most |M|^2 in the infinity norm. An error the Joseph
// step makes is carried by the later steps; one the prediction makes is carried
// first by its own step's Joseph factor, so the larger of the two is charged.
inline double covariance_premise(const analysis_maps &joseph, const analysis_maps &step)
{
    analysis_matrix later = analysis_matrix::Identity();
    double sum            = 0.0;
    for(std::size_t k = kGoldenSteps; k-- > 0;)
    {
        const double carried   = infinity_norm(later);
        const double predicted = infinity_norm(later * joseph[k]);
        sum += std::max(carried * carried, predicted * predicted);
        later = later * step[k];
    }
    return sum;
}

inline double covariance_surrogate(const analysis_matrix &posterior, const analysis_matrix &prior, const analysis_matrix &joseph, const Eigen::Vector2d &gain)
{
    const auto system           = estimation_demo<double>::plant();
    const auto settings         = estimation_demo<double>::config();
    const analysis_matrix a     = system.A.cwiseAbs();
    const analysis_matrix j     = joseph.cwiseAbs();
    const analysis_matrix ahead = a * posterior.cwiseAbs() * a.transpose() + settings.Q.cwiseAbs();
    const analysis_matrix after = j * prior.cwiseAbs() * j.transpose() + gain.cwiseAbs() * std::abs(settings.R(0, 0)) * gain.cwiseAbs().transpose();
    return std::max(ahead.maxCoeff(), after.maxCoeff());
}

inline double biquad_premise(const biquad_coeffs<double> &c)
{
    analysis_matrix f;
    f << -c.a1, 1.0, -c.a2, 0.0;
    analysis_matrix power = analysis_matrix::Identity();
    double sum            = 0.0;
    for(std::size_t m = 0; m < kGoldenSteps; ++m)
    {
        sum += infinity_norm(power);
        power = f * power;
    }
    return sum;
}

}

// The scale is the largest entry the covariance recursion forms in absolute
// values over the host run, each step's prior and gain re-formed from the
// posterior before it.
inline family_analysis analyze_estimation(estimation_demo<double> demo)
{
    const auto system   = estimation_demo<double>::plant();
    const auto settings = estimation_demo<double>::config();
    detail::analysis_maps joseph{};
    detail::analysis_maps step{};
    detail::analysis_matrix posterior = settings.P0;
    double scale                      = 0.0;
    double reference                  = 0.0;
    for(std::size_t k = 0; k < kGoldenSteps; ++k)
    {
        const detail::analysis_matrix prior = system.A * posterior * system.A.transpose() + settings.Q;
        const Eigen::Vector2d gain          = prior * system.C.transpose() / (system.C * prior * system.C.transpose() + settings.R)(0, 0);
        joseph[k]                           = detail::analysis_matrix::Identity() - gain * system.C;
        step[k]                             = joseph[k] * system.A;
        scale                               = std::max(scale, detail::covariance_surrogate(posterior, prior, joseph[k], gain));
        reference                           = demo.step();
        posterior                           = demo.covariance();
    }
    return {reference, scale, "largest covariance surrogate", scale, "propagated error norm sum", detail::covariance_premise(joseph, step), static_cast<double>(kGoldenSteps)};
}

inline family_analysis analyze_dsp(dsp_demo<double> demo)
{
    double peak      = 0.0;
    double reference = 0.0;
    for(std::int32_t k = 0; k < static_cast<std::int32_t>(kGoldenSteps); ++k)
    {
        reference = demo.step();
        peak      = std::max({peak, std::abs(dsp_demo<double>::input(k)), std::abs(reference)});
    }
    return {reference, biquad_scale(peak), "peak", peak, "power norm sum", detail::biquad_premise(demo.coefficients()), static_cast<double>(kGoldenSteps)};
}

inline family_analysis analyze_trajectory(trajectory_demo<double> demo)
{
    double reference = 0.0;
    for(std::size_t k = 0; k < kGoldenSteps; ++k)
        reference = demo.step();

    const auto command    = trajectory_demo<double>::command();
    const double duration = demo.profile().duration();
    const double scale    = double_s_scale(command.q0, command.q1, command.v_max, command.a_max, command.j_max, duration);
    return {reference, scale, "duration", duration, "limit not reached", demo.profile().is_degenerate() ? 1.0 : 0.0, 0.0};
}

}

#endif
