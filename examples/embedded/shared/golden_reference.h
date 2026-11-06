#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_GOLDEN_REFERENCE_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_GOLDEN_REFERENCE_H

#include "family_bounds.h"
#include "derived_tolerance.h"

#include <cstddef>
#include <cstdint>

namespace ctrlpp
{

namespace golden
{

constexpr double dt = 0.02;
constexpr std::int32_t steps = 200;

// The reference is a host double-precision build of the same algorithm (generate_golden.cpp), not an
// independent implementation: it shows a board's arithmetic agrees with a desktop build, and cannot catch a
// wrong Riccati solve, which the cross-validation suite under validation/ checks independently.
// The constants carry 17 significant digits, so each is the host's double exactly and the verdict charges
// the host run nothing for its gain.
constexpr double host_k0 = 9.4670936721545083;
constexpr double host_k1 = 5.281739638040154;
constexpr double host_final_norm = 2.0206867258019932e-05;

// The kernel runs one step per index from 0 through steps inclusive.
constexpr std::size_t sample_count = static_cast<std::size_t>(steps) + 1;
constexpr std::size_t nx = 2;
constexpr std::size_t nu = 1;

// Each family below is compared at the value its last step returns, against a host run of the same kernel
// at double precision. Every scale is a magnitude of that host run, never of a departure, so no constant
// here is fitted to the error it judges.
constexpr std::size_t kalman_nx = 2;
constexpr std::size_t kalman_ny = 1;
constexpr std::size_t kalman_plant_roundings = 1;

constexpr double host_kalman_velocity_variance = 1.3104323528650401;
constexpr double host_kalman_scale = 1.5685167138038882;
constexpr double host_biquad_output = -0.86508002466156531;
constexpr double host_biquad_scale = 16.492588552695299;
constexpr double host_trajectory_position = 6.744791666666667;
constexpr double host_trajectory_scale = 1206.09375;

// The predictive family runs the configuration the strict-zero allocation test pins, from x0 = (1, 0), for
// a warm-up and then the armed steps. Its reference is the exact finite-horizon law on the represented plant
// and weights.
constexpr std::size_t predictive_nx = 2;
constexpr std::size_t predictive_nu = 1;
constexpr std::size_t predictive_nh = 5;
constexpr int predictive_nv = static_cast<int>((predictive_nh + 1) * predictive_nx + predictive_nh * predictive_nu);
constexpr int predictive_max_m = static_cast<int>(predictive_nx * (predictive_nh + 1));
constexpr std::int32_t predictive_warmup_solves = 20;
constexpr std::size_t predictive_run_steps = static_cast<std::size_t>(predictive_warmup_solves) + sample_count;
constexpr double predictive_dt = 0.1;
constexpr double predictive_state_weight = 10.0;
constexpr double predictive_input_weight = 0.1;

constexpr double host_predictive_k0 = 1.9951805744562641;
constexpr double host_predictive_k1 = 6.4867948990510831;
constexpr double host_predictive_input_row = 33.820912149271386;
constexpr double host_predictive_kkt_condition = 2682.0131439072306;
constexpr double host_predictive_jacobian_sigma = 0.22307799946689291;
constexpr double host_predictive_jacobian_norm = 2.1000000000000001;
constexpr double host_predictive_cost = 189.50222139218934;

}

constexpr double estimation_roundings()
{
    const double per_step = covariance_roundings_per_filter_step(golden::kalman_nx, golden::kalman_ny, static_cast<double>(golden::kalman_plant_roundings));
    return static_cast<double>(golden::sample_count) * per_step;
}

constexpr double dsp_roundings()
{
    return static_cast<double>(golden::sample_count) * biquad_roundings_per_step();
}

// A sample's position depends on no earlier sample, so the run's length does not enter.
constexpr double trajectory_roundings()
{
    return double_s_roundings_per_sample();
}

template <class Scalar>
constexpr double estimation_tolerance()
{
    return paired_departure_bound<Scalar>(estimation_roundings(), golden::host_kalman_scale);
}

template <class Scalar>
constexpr double dsp_tolerance()
{
    return paired_departure_bound<Scalar>(dsp_roundings(), golden::host_biquad_scale);
}

template <class Scalar>
constexpr double trajectory_tolerance()
{
    return paired_departure_bound<Scalar>(trajectory_roundings(), golden::host_trajectory_scale);
}

}

#endif
