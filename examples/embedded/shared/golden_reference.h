#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_GOLDEN_REFERENCE_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_GOLDEN_REFERENCE_H

#include "family_bounds.h"
#include "derived_tolerance.h"

#include <cstddef>
#include <cstdint>

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

// The predictive family: the configuration the strict-zero allocation test
// pins, run from x0 = (1, 0) for a warm-up and then the armed steps. Its
// reference is the exact finite-horizon law on the represented plant and
// weights; the host constants below are that law's gain, the first-input row of
// the inverse KKT matrix with that matrix's condition number, the smallest
// singular value and the 2-norm bound of the constraint Jacobian, and the
// reference loop's realized cost.
constexpr std::size_t kPredictiveNx            = 2;
constexpr std::size_t kPredictiveNu            = 1;
constexpr std::size_t kPredictiveNh            = 5;
constexpr int kPredictiveNv                    = static_cast<int>((kPredictiveNh + 1) * kPredictiveNx + kPredictiveNh * kPredictiveNu);
constexpr int kPredictiveMaxM                  = static_cast<int>(kPredictiveNx * (kPredictiveNh + 1));
constexpr std::int32_t kPredictiveWarmupSolves = 20;
constexpr std::size_t kPredictiveRunSteps      = static_cast<std::size_t>(kPredictiveWarmupSolves) + kGoldenSteps;
constexpr double kPredictiveDt                 = 0.1;
constexpr double kPredictiveStateWeight        = 10.0;
constexpr double kPredictiveInputWeight        = 0.1;

constexpr double kHostPredictiveK0            = 1.9951805744562641;
constexpr double kHostPredictiveK1            = 6.4867948990510831;
constexpr double kHostPredictiveInputRow      = 33.820912149271386;
constexpr double kHostPredictiveKktCondition  = 2682.0131439072306;
constexpr double kHostPredictiveJacobianSigma = 0.22307799946689291;
constexpr double kHostPredictiveJacobianNorm  = 2.1000000000000001;
constexpr double kHostPredictiveCost          = 189.50222139218934;

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
