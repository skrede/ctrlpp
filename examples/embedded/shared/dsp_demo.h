#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_DSP_DEMO_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_DSP_DEMO_H

#include "golden_reference.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/dsp/biquad.h"
#include "ctrlpp/dsp/dsp_types.h"

#include <cstdint>
#include <numbers>

namespace ctrlpp {

inline const char *describe(dsp_error e)
{
    switch(e)
    {
        case dsp_error::non_positive_sample_rate:
            return "sample rate is not positive";
        case dsp_error::cutoff_exceeds_nyquist:
            return "cutoff is not below half the sample rate";
        case dsp_error::non_positive_q:
            return "quality factor is not positive";
        case dsp_error::non_positive_ripple:
            return "ripple is not positive";
        case dsp_error::non_finite_input:
            return "non-finite input";
    }
    return "unknown";
}

// biquad from its low-pass factory is what the embedded compile witness builds
// for this family, and the dsp row of the allocation matrix makes its claim
// about it.
template<class Scalar>
class dsp_demo
{
public:
    static ctrlpp::expected<dsp_demo, dsp_error> make()
    {
        const auto filter = biquad<Scalar>::low_pass(kCutoffHz, kSampleHz);
        if(!filter.has_value())
            return ctrlpp::unexpected(filter.error());
        return dsp_demo(*filter);
    }

    // A square wave below the cutoff, exact at either precision; the run ends on
    // the first sample after an edge, where the output is furthest from settled.
    static Scalar input(std::int32_t k)
    {
        return (k / kHalfPeriod) % 2 == 0 ? Scalar{1} : Scalar{-1};
    }

    Scalar step()
    {
        const Scalar x = input(index_);
        ++index_;
        return filter_.process(x);
    }

    const biquad_coeffs<Scalar> &coefficients() const
    {
        return filter_.coefficients();
    }

private:
    static constexpr Scalar kCutoffHz         = Scalar{100};
    static constexpr Scalar kSampleHz         = Scalar{1000};
    static constexpr std::int32_t kHalfPeriod = 25;

    // biquad_roundings_per_step charges the design's cosine and sine for a
    // normalized cutoff of at most one radian.
    static_assert(Scalar{2} * std::numbers::pi_v<Scalar> * kCutoffHz <= kSampleHz);

    biquad<Scalar> filter_;
    std::int32_t index_;

    explicit dsp_demo(biquad<Scalar> filter)
            : filter_(filter)
            , index_(0)
    {
    }
};

}

#endif
