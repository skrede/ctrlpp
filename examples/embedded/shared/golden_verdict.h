#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_GOLDEN_VERDICT_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_GOLDEN_VERDICT_H

#include "golden_reference.h"
#include "closed_loop_bound.h"
#include "control_loop_demo.h"
#include "derived_tolerance.h"

#include <array>
#include <cmath>
#include <limits>
#include <cstddef>

namespace ctrlpp {

using golden_bound = closed_loop_bound<2, 1, kGoldenSteps>;

// Every entry of A and B is formed from kDt in at most three roundings: B's
// first entry is dt*dt/2 from a rounded dt, and the halving is exact.
constexpr double kPlantRoundings = 3.0;

struct golden_verdict
{
    std::array<double, 2> gain_departure;
    double gain_tolerance;
    double norm_departure;
    double norm_tolerance;
    bool gain_pass;
    bool norm_pass;
    bool pass;
};

// Strict, so a departure exactly at its bound is refused.
inline bool within_bound(double departure, double bound)
{
    return departure < bound;
}

inline golden_bound make_golden_bound()
{
    control_loop_demo<double> reference = control_loop_demo<double>::plant();
    reference.K << kHostK0, kHostK1;
    return golden_bound(reference.A, reference.B, reference.K, reference.x);
}

// The board's gain and the host's are each within one counted bound of the
// exact design, so they may part by the sum of the two.
template<class Scalar>
double golden_gain_tolerance()
{
    const double gain_norm = std::abs(kHostK0) + std::abs(kHostK1);
    return gain_departure_bound<Scalar>(kGoldenNx, kGoldenNu, gain_norm) + gain_departure_bound<double>(kGoldenNx, kGoldenNu, gain_norm);
}

// The host run is the reference loop itself, so it carries only its own
// roundings. The 2-norm departure is at most sqrt(nx) times the infinity-norm
// one, and each side's norm evaluation adds nx + 1 roundings, charged at the
// board's own result.
template<class Scalar>
double golden_norm_tolerance(golden_bound &bound, double gain_departure, double final_norm)
{
    const double eps       = static_cast<double>(std::numeric_limits<Scalar>::epsilon());
    const double eps_host  = std::numeric_limits<double>::epsilon();
    const double board     = bound.final_state_departure(eps, kPlantRoundings, gain_departure);
    const double host      = bound.final_state_departure(eps_host, 0.0, 0.0);
    const double roundings = static_cast<double>(kGoldenNx + 1);
    const double scale     = std::sqrt(static_cast<double>(kGoldenNx));
    return scale * (board + host) + roundings * (eps * final_norm + eps_host * kHostFinalNorm);
}

// The gain departure feeding the trajectory bound is the board's own,
// measured, and admitted only once the gain check has passed against its
// counted bound. Each subtraction rounds once in double, charged by the factor.
template<class Scalar>
golden_verdict judge_golden(const control_loop_demo<Scalar> &demo, golden_bound &bound)
{
    golden_verdict verdict{};
    verdict.gain_departure[0] = static_cast<double>(demo.K(0, 0)) - kHostK0;
    verdict.gain_departure[1] = static_cast<double>(demo.K(0, 1)) - kHostK1;
    verdict.gain_tolerance    = golden_gain_tolerance<Scalar>();
    const double dev0         = std::abs(verdict.gain_departure[0]);
    const double dev1         = std::abs(verdict.gain_departure[1]);
    verdict.gain_pass         = within_bound(dev0, verdict.gain_tolerance) && within_bound(dev1, verdict.gain_tolerance);

    const double measured   = (dev0 + dev1) * (1.0 + std::numeric_limits<double>::epsilon());
    const double final_norm = static_cast<double>(demo.x.norm());
    verdict.norm_departure  = std::abs(final_norm - kHostFinalNorm);
    verdict.norm_tolerance  = golden_norm_tolerance<Scalar>(bound, measured, final_norm);
    verdict.norm_pass       = within_bound(verdict.norm_departure, verdict.norm_tolerance);
    verdict.pass            = verdict.gain_pass && verdict.norm_pass;
    return verdict;
}

}

#endif
