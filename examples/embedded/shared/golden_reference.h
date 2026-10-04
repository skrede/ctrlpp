#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_GOLDEN_REFERENCE_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_GOLDEN_REFERENCE_H

#include "family_bounds.h"
#include "derived_tolerance.h"

#include <cstddef>

namespace ctrlpp
{

constexpr double kDt    = 0.02;
constexpr int    kSteps = 200;

// The reference below is an independent host BUILD of the same algorithm at
// double precision (generate_golden.cpp), not an independent implementation. It
// establishes that a board's arithmetic and its FPU agree with a desktop double
// build on this problem, which is a portability claim; it cannot detect a wrong
// Riccati solve, because both sides would be wrong identically. The
// independent-implementation check is the cross-validation suite under
// validation/. The constants carry 17 significant digits so each is the host's
// double exactly, which is what lets the verdict charge the host run nothing
// for its gain.
constexpr double kHostK0        = 9.4670936721545083;
constexpr double kHostK1        = 5.281739638040154;
constexpr double kHostFinalNorm = 2.0206867258019932e-05;

// The kernel runs one step per index from 0 through kSteps inclusive.
constexpr std::size_t kGoldenSteps = static_cast<std::size_t>(kSteps) + 1;
constexpr std::size_t kGoldenNx    = 2;
constexpr std::size_t kGoldenNu    = 1;

// Each family below is compared at the value its last step returns, against a
// host run of the same kernel at double precision. Every scale is a magnitude
// of that host run, never of a departure, so no constant here is fitted to the
// error it judges.
constexpr std::size_t kKalmanNx             = 2;
constexpr std::size_t kKalmanNy             = 1;
constexpr std::size_t kKalmanPlantRoundings = 1;

constexpr double kHostKalmanVelocityVariance = 1.3104323528650401;
constexpr double kHostKalmanScale            = 1.5685167138038882;
constexpr double kHostBiquadOutput           = -0.86508002466156553;
constexpr double kHostBiquadScale            = 16.492588552695306;
constexpr double kHostTrajectoryPosition     = 6.744791666666667;
constexpr double kHostTrajectoryScale        = 1206.09375;

constexpr double estimation_roundings()
{
    const double per_step = covariance_roundings_per_filter_step(kKalmanNx, kKalmanNy, static_cast<double>(kKalmanPlantRoundings));
    return static_cast<double>(kGoldenSteps) * per_step;
}

constexpr double dsp_roundings()
{
    return static_cast<double>(kGoldenSteps) * biquad_roundings_per_step();
}

// A sample's position depends on no earlier sample, so the run's length does
// not enter.
constexpr double trajectory_roundings()
{
    return double_s_roundings_per_sample();
}

template<class Scalar>
constexpr double estimation_tolerance()
{
    return paired_departure_bound<Scalar>(estimation_roundings(), kHostKalmanScale);
}

template<class Scalar>
constexpr double dsp_tolerance()
{
    return paired_departure_bound<Scalar>(dsp_roundings(), kHostBiquadScale);
}

template<class Scalar>
constexpr double trajectory_tolerance()
{
    return paired_departure_bound<Scalar>(trajectory_roundings(), kHostTrajectoryScale);
}

}

#endif
