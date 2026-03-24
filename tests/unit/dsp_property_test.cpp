#include "ctrlpp/dsp/biquad.h"
#include "ctrlpp/dsp/fir.h"

#include <catch2/catch_test_macros.hpp>
#include <rapidcheck.h>
#include <rapidcheck/catch.h>

#include <array>
#include <cmath>
#include <cstddef>

using namespace ctrlpp;

namespace
{

auto bounded_double(double lo, double hi) -> rc::Gen<double>
{
    return rc::gen::map(rc::gen::inRange(0, 1000000), [lo, hi](int x) { return lo + (hi - lo) * (static_cast<double>(x) / 1000000.0); });
}

auto gen_signal(std::size_t len) -> rc::Gen<std::vector<double>>
{
    return rc::gen::container<std::vector<double>>(len, bounded_double(-10.0, 10.0));
}

} // namespace

TEST_CASE("dsp property tests", "[dsp][property]")
{
    SECTION("biquad low-pass is linear: filter(a*x + b*y) == a*filter(x) + b*filter(y)")
    {
        rc::prop("biquad linearity", [](void)
                 {
            constexpr std::size_t N = 50;
            auto x_sig = *gen_signal(N);
            auto y_sig = *gen_signal(N);
            auto a = *bounded_double(-3.0, 3.0);
            auto b = *bounded_double(-3.0, 3.0);

            auto f1 = biquad<double>::low_pass(10.0, 100.0);
            auto f2 = biquad<double>::low_pass(10.0, 100.0);
            auto f3 = biquad<double>::low_pass(10.0, 100.0);

            double max_err = 0.0;
            for(std::size_t i = 0; i < N; ++i)
            {
                double y_ax = f1.process(a * x_sig[i] + b * y_sig[i]);
                double y_a = f2.process(x_sig[i]);
                double y_b = f3.process(y_sig[i]);
                double combined = a * y_a + b * y_b;
                max_err = std::max(max_err, std::abs(y_ax - combined));
            }

            RC_ASSERT(max_err < 1e-10); });
    }

    SECTION("fir is linear: fir(a*x + b*y) == a*fir(x) + b*fir(y)")
    {
        rc::prop("fir linearity", [](void)
                 {
            constexpr std::size_t N = 50;
            constexpr std::size_t NTAPS = 5;
            auto x_sig = *gen_signal(N);
            auto y_sig = *gen_signal(N);
            auto a = *bounded_double(-3.0, 3.0);
            auto b = *bounded_double(-3.0, 3.0);

            std::array<double, NTAPS> taps{0.2, 0.2, 0.2, 0.2, 0.2};
            fir<double, NTAPS> f1(taps);
            fir<double, NTAPS> f2(taps);
            fir<double, NTAPS> f3(taps);

            double max_err = 0.0;
            for(std::size_t i = 0; i < N; ++i)
            {
                double y_ax = f1.process(a * x_sig[i] + b * y_sig[i]);
                double y_a = f2.process(x_sig[i]);
                double y_b = f3.process(y_sig[i]);
                double combined = a * y_a + b * y_b;
                max_err = std::max(max_err, std::abs(y_ax - combined));
            }

            RC_ASSERT(max_err < 1e-10); });
    }

    SECTION("biquad low-pass DC gain is unity")
    {
        rc::prop("DC gain unity", [](void)
                 {
            auto dc_val = *bounded_double(0.1, 10.0);

            auto filter = biquad<double>::low_pass(20.0, 1000.0);

            double output = 0.0;
            for(int i = 0; i < 500; ++i)
                output = filter.process(dc_val);

            // After settling, output should approach dc_val
            RC_ASSERT(std::abs(output - dc_val) < 0.01 * dc_val + 1e-10); });
    }
}
