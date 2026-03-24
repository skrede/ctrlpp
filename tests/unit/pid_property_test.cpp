#include "ctrlpp/control/pid.h"
#include "ctrlpp/control/pid_policies.h"

#include <catch2/catch_test_macros.hpp>
#include <rapidcheck/catch.h>
#include <rapidcheck.h>

#include <cmath>
#include <cstddef>

namespace
{

using Vec1 = ctrlpp::Vector<double, 1>;

auto vec1(double v) -> Vec1
{
    Vec1 r;
    r << v;
    return r;
}

auto bounded_double(double lo, double hi) -> rc::Gen<double>
{
    return rc::gen::map(rc::gen::inRange(0, 1000000), [lo, hi](int x) { return lo + (hi - lo) * (static_cast<double>(x) / 1000000.0); });
}

} // namespace

TEST_CASE("pid bibo stability property", "[pid][property]")
{
    rc::prop("output is always within clamped bounds",
             []
             {
                 auto kp = *bounded_double(0.01, 100.0);
                 auto ki = *bounded_double(0.01, 100.0);
                 auto kd = *bounded_double(0.01, 100.0);
                 auto out_min = *bounded_double(-200.0, -1.0);
                 auto out_max = *bounded_double(1.0, 200.0);
                 auto n_steps = *rc::gen::inRange(10, 200);

                 using PidType = ctrlpp::pid<double, 1, 1, 1, ctrlpp::anti_windup<ctrlpp::clamping>>;
                 typename PidType::config_type cfg{};
                 cfg.kp = vec1(kp);
                 cfg.ki = vec1(ki);
                 cfg.kd = vec1(kd);
                 cfg.output_min = vec1(out_min);
                 cfg.output_max = vec1(out_max);

                 PidType pid(cfg);

                 constexpr double dt = 0.01;
                 constexpr double eps = 1e-9;

                 for(int i = 0; i < n_steps; ++i)
                 {
                     auto sp = *bounded_double(-100.0, 100.0);
                     auto pv = *bounded_double(-100.0, 100.0);
                     auto u = pid.compute(vec1(sp), vec1(pv), dt);
                     RC_ASSERT(u[0] >= out_min - eps);
                     RC_ASSERT(u[0] <= out_max + eps);
                 }
             });
}

TEST_CASE("pid proportional-only constant error", "[pid][property]")
{
    rc::prop("P-only output equals kp times error",
             []
             {
                 auto kp = *bounded_double(0.01, 100.0);

                 using PidType = ctrlpp::pid<double, 1, 1, 1>;
                 typename PidType::config_type cfg{};
                 cfg.kp = vec1(kp);

                 PidType pid(cfg);

                 auto target = *bounded_double(-100.0, 100.0);
                 constexpr double dt = 0.01;
                 constexpr double tol = 1e-9;

                 // First step: error = target - 0 = target, output = kp * target
                 auto u = pid.compute(vec1(target), vec1(0.0), dt);
                 RC_ASSERT(std::abs(u[0] - kp * target) < tol);

                 // Subsequent steps with same input: P contribution unchanged,
                 // ki=kd=0 so no integral or derivative
                 for(int i = 0; i < 10; ++i)
                 {
                     u = pid.compute(vec1(target), vec1(0.0), dt);
                     RC_ASSERT(std::abs(u[0] - kp * target) < tol);
                 }
             });
}

TEST_CASE("pid robustness - extreme inputs do not crash", "[pid][property]")
{
    rc::prop("extreme inputs produce finite output",
             []
             {
                 auto kp = *bounded_double(0.0, 1e6);
                 auto ki = *bounded_double(0.0, 1e6);
                 auto kd = *bounded_double(0.0, 1e6);
                 auto dt = *bounded_double(1e-10, 100.0);
                 auto sp = *bounded_double(-1e6, 1e6);
                 auto pv = *bounded_double(-1e6, 1e6);

                 using PidType = ctrlpp::pid<double, 1, 1, 1, ctrlpp::anti_windup<ctrlpp::clamping>>;
                 typename PidType::config_type cfg{};
                 cfg.kp = vec1(kp);
                 cfg.ki = vec1(ki);
                 cfg.kd = vec1(kd);
                 cfg.output_min = vec1(-1e9);
                 cfg.output_max = vec1(1e9);

                 PidType pid(cfg);

                 auto n_steps = *rc::gen::inRange(1, 50);
                 for(int i = 0; i < n_steps; ++i)
                 {
                     auto sp_i = *bounded_double(-1e6, 1e6);
                     auto pv_i = *bounded_double(-1e6, 1e6);
                     auto u = pid.compute(vec1(sp_i), vec1(pv_i), dt);
                     RC_ASSERT(std::isfinite(u[0]));
                 }
             });
}
