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

namespace ctrlpp
{

using predictive_workspace = closed_loop_bound<static_cast<int>(golden::predictive_nx), static_cast<int>(golden::predictive_nu), golden::predictive_run_steps>;

struct predictive_record
{
    std::array<double, golden::predictive_run_steps> certificate;
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
    const double dt = golden::predictive_dt;
    Eigen::Matrix2d A;
    A << 1.0, dt, 0.0, 1.0;
    const Eigen::Vector2d B(0.0, dt);
    Eigen::RowVector2d K;
    K << golden::host_predictive_k0, golden::host_predictive_k1;
    return predictive_workspace(A, B, K, Eigen::Vector2d(1.0, 0.0));
}

// The realized closed-loop cost, the quantity the controller is built to keep
// small; every term is non-negative, so its roundings are charged at the sum.
inline double stage_cost(const Eigen::Vector2d& x, double u)
{
    return golden::predictive_state_weight * x.squaredNorm() + golden::predictive_input_weight * (u * u);
}

constexpr double cost_roundings()
{
    return 2.0 * static_cast<double>(golden::predictive_nx) + 3.0 + static_cast<double>(golden::predictive_run_steps);
}

// Runs between solves, so nothing here falls inside an armed window. The gate
// admits the departure the certificate allows from the exact law, the host
// gain's own departure from it, and the four roundings of the comparison.
inline void observe(predictive_record& record, const Eigen::Vector2d& before, const predictive_demo<double>& demo, std::size_t k)
{
    const double eps = std::numeric_limits<double>::epsilon();
    const double u = demo.input()[0];
    const double x_inf = before.lpNorm<Eigen::Infinity>();
    const double z_inf = demo.controller().last_solution().lpNorm<Eigen::Infinity>();
    const double certified = input_certificate(z_inf, argmin_settings<double>{}.kkt_tol);
    const double gain = std::abs(golden::host_predictive_k0) + std::abs(golden::host_predictive_k1);
    const double law = static_cast<double>(golden::predictive_nx) * predictive_gain_departure() * x_inf;
    const double gate = certified + law + counted_departure_bound<double>(4.0, std::abs(u) + gain * x_inf);
    const double departure = std::abs(u + golden::host_predictive_k0 * before[0] + golden::host_predictive_k1 * before[1]);
    record.certificate[k] = certified * (1.0 + 2.0 * eps);
    record.worst_gate = std::max(record.worst_gate, departure / gate);
    record.cost += stage_cost(before, u);
    if(demo.controller().diagnostics().stop_criterion == nlp_stop_criterion::stationarity)
        ++record.stationarity_stops;
}

// The first golden::predictive_warmup_solves steps walk the solver into steady state;
// the rest are the window the board's figures cover.
constexpr bool in_predictive_window(std::size_t k)
{
    return k >= static_cast<std::size_t>(golden::predictive_warmup_solves);
}

// The hooks bracket the solve and the plant step alone, on every step, and are
// told which step it is. A refused construction returns before any hook runs,
// since the hooks may read the controller.
template <class Enter, class Leave>
ctrlpp::expected<void, solver_error> run_predictive(predictive_demo<double>& demo, predictive_record& record, Enter enter, Leave leave)
{
    record = predictive_record{{}, 0.0, 0.0, 0};
    if(!demo.constructed())
        return ctrlpp::unexpected(solver_error::setup_incomplete);
    for(std::size_t k = 0; k < golden::predictive_run_steps; ++k)
    {
        const Eigen::Vector2d before = demo.state();
        enter(k);
        const auto advanced = demo.advance();
        leave(k);
        if(!advanced.has_value())
            return advanced;
        observe(record, before, demo, k);
    }
    return {};
}

// The departure of a run's cost from the reference loop's, for state departures
// d and input departures v at reference state norms X: with Q = q*I,
// |x'Qx - y'Qy| <= nx*q*d*(2*X + d), and likewise r*v*(2*|K|*X + v) for R.
template <class InputDeparture>
double cost_departure(const predictive_workspace& workspace, const std::array<double, golden::predictive_run_steps + 1>& d, InputDeparture input)
{
    const auto& X = workspace.reference_state_norms();
    const double nx = static_cast<double>(golden::predictive_nx);
    const double gain = workspace.gain_norm();
    double total = 0.0;
    for(std::size_t k = 0; k < golden::predictive_run_steps; ++k)
    {
        const double v = input(k, X[k], d[k]);
        total += nx * golden::predictive_state_weight * d[k] * (2.0 * X[k] + d[k]) + golden::predictive_input_weight * v * (2.0 * gain * X[k] + v);
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
inline double host_cost_departure(predictive_workspace& workspace)
{
    const double eps = std::numeric_limits<double>::epsilon();
    const double gain = workspace.gain_norm();
    const double roundings = static_cast<double>(2 * golden::predictive_nx - 1);
    constexpr std::array<double, golden::predictive_run_steps> none{};
    const auto input = [&](std::size_t, double X, double d) { return gain * d + counted_departure_bound<double>(roundings, gain * (X + d)); };
    return cost_departure(workspace, workspace.input_departures(eps, none, 0.0), input) + cost_rounding(golden::host_predictive_cost);
}

// The board's run is charged the certified departure of every solve from the
// exact law, propagated through the reference loop; the certificate is a bound
// only under the premise, and the gate is where it is seen to hold.
inline predictive_verdict judge_predictive(const predictive_record& record, predictive_workspace& workspace)
{
    const double eps = std::numeric_limits<double>::epsilon();
    const double gain = workspace.gain_norm();
    const double law = static_cast<double>(golden::predictive_nx) * predictive_gain_departure() * (1.0 + 2.0 * eps);
    const double host = host_cost_departure(workspace);
    const auto board_input = [&](std::size_t k, double X, double d) { return gain * d + record.certificate[k] + law * (X + d); };
    const double board = cost_departure(workspace, workspace.input_departures(eps, record.certificate, law), board_input);

    predictive_verdict verdict{};
    verdict.premise = record.stationarity_stops == static_cast<std::uint32_t>(golden::predictive_run_steps);
    verdict.worst_gate = record.worst_gate;
    verdict.gate = record.worst_gate < 1.0;
    verdict.cost = record.cost;
    verdict.departure = std::abs(record.cost - golden::host_predictive_cost);
    verdict.bound = board + host + cost_rounding(record.cost);
    verdict.pass = verdict.premise && verdict.gate && within_bound(verdict.departure, verdict.bound);
    return verdict;
}

}

#endif
