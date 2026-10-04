#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_PREDICTIVE_VERDICT_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_PREDICTIVE_VERDICT_H

#include "golden_verdict.h"
#include "predictive_demo.h"
#include "golden_reference.h"
#include "predictive_bound.h"
#include "closed_loop_bound.h"
#include "derived_tolerance.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/mpc/nlp_types.h"
#include "ctrlpp/mpc/argmin_policies.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <algorithm>

namespace ctrlpp {

using predictive_workspace = closed_loop_bound<static_cast<int>(kPredictiveNx), static_cast<int>(kPredictiveNu), kPredictiveRunSteps>;

struct predictive_record
{
    std::array<double, kPredictiveRunSteps> certificate;
    double cost;
    double worst_gate;
    std::uint32_t stationarity_stops;
};

struct predictive_verdict
{
    bool premise;
    bool gate;
    double worst_gate;
    double cost;
    double departure;
    double bound;
    bool pass;
};

inline predictive_workspace make_predictive_workspace()
{
    const double dt = kPredictiveDt;
    Eigen::Matrix2d A;
    A << 1.0, dt, 0.0, 1.0;
    const Eigen::Vector2d B(0.0, dt);
    Eigen::RowVector2d K;
    K << kHostPredictiveK0, kHostPredictiveK1;
    return predictive_workspace(A, B, K, Eigen::Vector2d(1.0, 0.0));
}

// The realized closed-loop cost, the quantity the controller is built to keep
// small; every term is non-negative, so its roundings are charged at the sum.
inline double stage_cost(const Eigen::Vector2d &x, double u)
{
    return kPredictiveStateWeight * x.squaredNorm() + kPredictiveInputWeight * (u * u);
}

constexpr double cost_roundings()
{
    return 2.0 * static_cast<double>(kPredictiveNx) + 3.0 + static_cast<double>(kPredictiveRunSteps);
}

// Runs between solves, so nothing here falls inside an armed window. The gate
// admits the departure the certificate allows from the exact law, the host
// gain's own departure from it, and the four roundings of the comparison.
inline void observe(predictive_record &record, const Eigen::Vector2d &before, const predictive_demo<double> &demo, std::size_t k)
{
    const double eps       = std::numeric_limits<double>::epsilon();
    const double u         = demo.input()[0];
    const double x_inf     = before.lpNorm<Eigen::Infinity>();
    const double z_inf     = demo.controller().last_solution().lpNorm<Eigen::Infinity>();
    const double certified = input_certificate(z_inf, argmin_settings<double>{}.kkt_tol);
    const double gain      = std::abs(kHostPredictiveK0) + std::abs(kHostPredictiveK1);
    const double law       = static_cast<double>(kPredictiveNx) * predictive_gain_departure() * x_inf;
    const double gate      = certified + law + counted_departure_bound<double>(4.0, std::abs(u) + gain * x_inf);
    const double departure = std::abs(u + kHostPredictiveK0 * before[0] + kHostPredictiveK1 * before[1]);
    record.certificate[k]  = certified * (1.0 + 2.0 * eps);
    record.worst_gate      = std::max(record.worst_gate, departure / gate);
    record.cost += stage_cost(before, u);
    if(demo.controller().diagnostics().stop_criterion == nlp_stop_criterion::stationarity)
        ++record.stationarity_stops;
}

// The first kPredictiveWarmupSolves steps walk the solver into steady state
// unarmed; the rest are armed around the solve and the plant step alone.
template<class Arm, class Disarm>
ctrlpp::expected<void, solver_error> run_predictive(predictive_demo<double> &demo, predictive_record &record, Arm arm, Disarm disarm)
{
    record = predictive_record{{}, 0.0, 0.0, 0};
    for(std::size_t k = 0; k < kPredictiveRunSteps; ++k)
    {
        const bool armed             = k >= static_cast<std::size_t>(kPredictiveWarmupSolves);
        const Eigen::Vector2d before = demo.state();
        if(armed)
            arm();
        const auto advanced = demo.advance();
        if(armed)
            disarm();
        if(!advanced.has_value())
            return advanced;
        observe(record, before, demo, k);
    }
    return {};
}

// The departure of a run's cost from the reference loop's, for state departures
// d and input departures v at reference state norms X: with Q = q*I,
// |x'Qx - y'Qy| <= nx*q*d*(2*X + d), and likewise r*v*(2*|K|*X + v) for R.
template<class InputDeparture>
double cost_departure(const predictive_workspace &workspace, const std::array<double, kPredictiveRunSteps + 1> &d, InputDeparture input)
{
    const auto &X     = workspace.reference_state_norms();
    const double nx   = static_cast<double>(kPredictiveNx);
    const double gain = workspace.gain_norm();
    double total      = 0.0;
    for(std::size_t k = 0; k < kPredictiveRunSteps; ++k)
    {
        const double v = input(k, X[k], d[k]);
        total += nx * kPredictiveStateWeight * d[k] * (2.0 * X[k] + d[k]) + kPredictiveInputWeight * v * (2.0 * gain * X[k] + v);
    }
    return total;
}

inline double cost_rounding(double cost)
{
    const double eps = std::numeric_limits<double>::epsilon();
    return counted_departure_bound<double>(cost_roundings(), cost) / (1.0 - cost_roundings() * eps);
}

// The host's reference run carries only its own roundings: those of its
// states, the 2*nx - 1 of each u = -K*x, and those of its cost.
inline double host_cost_departure(predictive_workspace &workspace)
{
    const double eps       = std::numeric_limits<double>::epsilon();
    const double gain      = workspace.gain_norm();
    const double roundings = static_cast<double>(2 * kPredictiveNx - 1);
    constexpr std::array<double, kPredictiveRunSteps> none{};
    const auto input = [&](std::size_t, double X, double d) { return gain * d + counted_departure_bound<double>(roundings, gain * (X + d)); };
    return cost_departure(workspace, workspace.input_departures(eps, none, 0.0), input) + cost_rounding(kHostPredictiveCost);
}

// The board's run is charged the certified departure of every solve from the
// exact law, propagated through the reference loop; the certificate is a bound
// only under the premise, and the gate is where it is seen to hold.
inline predictive_verdict judge_predictive(const predictive_record &record, predictive_workspace &workspace)
{
    const double eps       = std::numeric_limits<double>::epsilon();
    const double gain      = workspace.gain_norm();
    const double law       = static_cast<double>(kPredictiveNx) * predictive_gain_departure() * (1.0 + 2.0 * eps);
    const double host      = host_cost_departure(workspace);
    const auto board_input = [&](std::size_t k, double X, double d) { return gain * d + record.certificate[k] + law * (X + d); };
    const double board     = cost_departure(workspace, workspace.input_departures(eps, record.certificate, law), board_input);

    predictive_verdict verdict{};
    verdict.premise    = record.stationarity_stops == static_cast<std::uint32_t>(kPredictiveRunSteps);
    verdict.worst_gate = record.worst_gate;
    verdict.gate       = record.worst_gate < 1.0;
    verdict.cost       = record.cost;
    verdict.departure  = std::abs(record.cost - kHostPredictiveCost);
    verdict.bound      = board + host + cost_rounding(record.cost);
    verdict.pass       = verdict.premise && verdict.gate && within_bound(verdict.departure, verdict.bound);
    return verdict;
}

}

#endif
