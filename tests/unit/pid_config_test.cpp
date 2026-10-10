#include "hardening_helpers.h"
#include "ctrlpp/control/pid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include <limits>
#include <vector>
#include <cstddef>
#include <cstring>
#include <type_traits>

namespace
{

using ctrlpp::pid_config_error;

template <typename Values>
auto same_bits(const Values& a, const Values& b) -> bool
{
    return std::memcmp(a.data(), b.data(), sizeof(a[0]) * static_cast<std::size_t>(a.size())) == 0;
}

template <typename Scalar>
auto not_a_number() -> Scalar
{
    return std::numeric_limits<Scalar>::quiet_NaN();
}

template <typename Scalar>
auto infinity() -> Scalar
{
    return std::numeric_limits<Scalar>::infinity();
}

template <typename Scalar>
using filtered_pid = ctrlpp::pid<Scalar, 2, ctrlpp::isa_form, ctrlpp::setpoint_filter, ctrlpp::pv_filter, ctrlpp::rate_limit,
                                 ctrlpp::anti_windup<ctrlpp::back_calc>, ctrlpp::deriv_filter>;

template <typename Scalar>
using conditional_pid = ctrlpp::pid<Scalar, 2, ctrlpp::isa_form, ctrlpp::rate_limit,
                                    ctrlpp::anti_windup<ctrlpp::conditional_integration>, ctrlpp::deriv_filter>;

template <typename Scalar>
using parallel_pid = ctrlpp::pid<Scalar, 2, ctrlpp::setpoint_filter, ctrlpp::pv_filter, ctrlpp::rate_limit,
                                 ctrlpp::anti_windup<ctrlpp::back_calc>, ctrlpp::deriv_filter>;

template <typename Config>
auto with_base_fields(Config cfg) -> Config
{
    using Scalar = typename Config::vector_t::Scalar;
    cfg.kp = Config::vector_t::Constant(Scalar{2});
    cfg.ki = Config::vector_t::Constant(static_cast<Scalar>(0.5));
    cfg.kd = Config::vector_t::Constant(static_cast<Scalar>(0.05));
    cfg.output_min = Config::vector_t::Constant(Scalar{-10});
    cfg.output_max = Config::vector_t::Constant(Scalar{10});
    return cfg;
}

template <typename Scalar>
auto filtered_baseline() -> typename filtered_pid<Scalar>::config_type
{
    auto cfg = with_base_fields(typename filtered_pid<Scalar>::config_type{});
    cfg.template policy<ctrlpp::setpoint_filter>().tf = {static_cast<Scalar>(0.05), static_cast<Scalar>(0.05)};
    cfg.template policy<ctrlpp::pv_filter>().tf = {static_cast<Scalar>(0.02), static_cast<Scalar>(0.02)};
    cfg.template policy<ctrlpp::rate_limit>().rate_max = {Scalar{50}, Scalar{50}};
    cfg.template policy<ctrlpp::deriv_filter>().n = {Scalar{10}, Scalar{10}};
    return cfg;
}

template <typename Scalar>
auto conditional_baseline() -> typename conditional_pid<Scalar>::config_type
{
    auto cfg = with_base_fields(typename conditional_pid<Scalar>::config_type{});
    cfg.template policy<ctrlpp::rate_limit>().rate_max = {Scalar{50}, Scalar{50}};
    cfg.template policy<ctrlpp::anti_windup<ctrlpp::conditional_integration>>().error_threshold = {Scalar{2}, Scalar{2}};
    cfg.template policy<ctrlpp::deriv_filter>().n = {Scalar{10}, Scalar{10}};
    return cfg;
}

template <typename Config>
struct fault
{
    void (*apply)(Config&);
    pid_config_error error;
};

template <typename Config>
void gain_fault(Config& cfg)
{
    cfg.kp[1] = not_a_number<typename Config::vector_t::Scalar>();
}

template <typename Config>
void zero_integral_time_fault(Config& cfg)
{
    cfg.ki[1] = typename Config::vector_t::Scalar{0};
}

template <typename Config>
void limits_fault(Config& cfg)
{
    using Scalar = typename Config::vector_t::Scalar;
    cfg.output_min[1] = Scalar{1};
    cfg.output_max[1] = Scalar{-1};
}

template <typename Config>
void weight_fault(Config& cfg)
{
    cfg.b[1] = infinity<typename Config::vector_t::Scalar>();
}

template <typename Config>
void unrepresentable_fault(Config& cfg)
{
    using Scalar = typename Config::vector_t::Scalar;
    cfg.kp[1] = std::numeric_limits<Scalar>::max();
    cfg.kd[1] = Scalar{2};
}

// The faults in check order, each confined to the second channel, so every check is also shown to
// scan past the first channel.
template <typename Scalar>
auto filtered_faults() -> std::vector<fault<typename filtered_pid<Scalar>::config_type>>
{
    using C = typename filtered_pid<Scalar>::config_type;
    return {
        {gain_fault<C>, pid_config_error::non_finite_gain},
        {zero_integral_time_fault<C>, pid_config_error::zero_integral_time},
        {limits_fault<C>, pid_config_error::invalid_output_limits},
        {weight_fault<C>, pid_config_error::non_finite_setpoint_weight},
        {[](C& c) { c.template policy<ctrlpp::pv_filter>().tf[1] = Scalar{-1}; }, pid_config_error::invalid_filter_time_constant},
        {[](C& c) { c.template policy<ctrlpp::rate_limit>().rate_max[1] = not_a_number<Scalar>(); }, pid_config_error::invalid_rate_limit},
        {[](C& c) { c.template policy<ctrlpp::anti_windup<ctrlpp::back_calc>>().kb[1] = infinity<Scalar>(); },
         pid_config_error::invalid_back_calculation_gain},
        {[](C& c) { c.template policy<ctrlpp::deriv_filter>().n[1] = not_a_number<Scalar>(); }, pid_config_error::invalid_derivative_filter},
        {[](C& c) { c.ki[1] = static_cast<Scalar>(-0.5); }, pid_config_error::undefined_tracking_time},
        {[](C& c) { c.template policy<ctrlpp::deriv_filter>().n[1] = Scalar{0}; }, pid_config_error::undefined_derivative_filter},
        {unrepresentable_fault<C>, pid_config_error::unrepresentable_gain},
    };
}

template <typename Scalar>
auto conditional_faults() -> std::vector<fault<typename conditional_pid<Scalar>::config_type>>
{
    using C = typename conditional_pid<Scalar>::config_type;
    return {
        {gain_fault<C>, pid_config_error::non_finite_gain},
        {zero_integral_time_fault<C>, pid_config_error::zero_integral_time},
        {limits_fault<C>, pid_config_error::invalid_output_limits},
        {weight_fault<C>, pid_config_error::non_finite_setpoint_weight},
        {[](C& c) { c.template policy<ctrlpp::rate_limit>().rate_max[1] = Scalar{-1}; }, pid_config_error::invalid_rate_limit},
        {[](C& c) { c.template policy<ctrlpp::anti_windup<ctrlpp::conditional_integration>>().error_threshold[1] = not_a_number<Scalar>(); },
         pid_config_error::invalid_error_threshold},
        {[](C& c) { c.template policy<ctrlpp::deriv_filter>().n[1] = infinity<Scalar>(); }, pid_config_error::invalid_derivative_filter},
        {[](C& c) { c.template policy<ctrlpp::deriv_filter>().n[1] = Scalar{-1}; }, pid_config_error::undefined_derivative_filter},
        {unrepresentable_fault<C>, pid_config_error::unrepresentable_gain},
    };
}

template <typename Pid>
void require_refusal_order(const typename Pid::config_type& baseline, const std::vector<fault<typename Pid::config_type>>& faults)
{
    REQUIRE(Pid::create(baseline).has_value());
    for(std::size_t i = 0; i < faults.size(); ++i)
    {
        for(std::size_t j = i; j < faults.size(); ++j)
        {
            CAPTURE(i, j);
            auto cfg = baseline;
            faults[j].apply(cfg);
            faults[i].apply(cfg);
            auto const created = Pid::create(cfg);
            REQUIRE_FALSE(created.has_value());
            REQUIRE(created.error() == faults[i].error);
        }
    }
}

template <typename Pid>
auto refusal_of(const typename Pid::config_type& cfg) -> pid_config_error
{
    auto const created = Pid::create(cfg);
    REQUIRE_FALSE(created.has_value());
    return created.error();
}

template <typename Pid>
auto drive(Pid& controller, int cycles) -> void
{
    using vector_t = typename Pid::vector_t;
    using Scalar = typename vector_t::Scalar;
    for(int k = 0; k < cycles; ++k)
    {
        auto const sp = vector_t::Constant(static_cast<Scalar>(1 + k % 3));
        auto const meas = vector_t::Constant(static_cast<Scalar>(k) * static_cast<Scalar>(0.25));
        REQUIRE(controller.compute(sp, meas, static_cast<Scalar>(0.01)).has_value());
    }
}

template <typename Config>
auto same_config(const Config& a, const Config& b) -> bool
{
    return same_bits(a.kp, b.kp) && same_bits(a.ki, b.ki) && same_bits(a.kd, b.kd) && same_bits(a.output_min, b.output_min)
        && same_bits(a.output_max, b.output_max) && same_bits(a.b, b.b) && same_bits(a.c, b.c)
        && a.derivative_on_error == b.derivative_on_error;
}

template <typename Scalar>
auto same_policies(const typename filtered_pid<Scalar>::config_type& a, const typename filtered_pid<Scalar>::config_type& b) -> bool
{
    return same_bits(a.template policy<ctrlpp::setpoint_filter>().tf, b.template policy<ctrlpp::setpoint_filter>().tf)
        && same_bits(a.template policy<ctrlpp::pv_filter>().tf, b.template policy<ctrlpp::pv_filter>().tf)
        && same_bits(a.template policy<ctrlpp::rate_limit>().rate_max, b.template policy<ctrlpp::rate_limit>().rate_max)
        && same_bits(a.template policy<ctrlpp::anti_windup<ctrlpp::back_calc>>().kb,
                     b.template policy<ctrlpp::anti_windup<ctrlpp::back_calc>>().kb)
        && same_bits(a.template policy<ctrlpp::deriv_filter>().n, b.template policy<ctrlpp::deriv_filter>().n)
        && a.template policy<ctrlpp::deriv_filter>().tf == b.template policy<ctrlpp::deriv_filter>().tf;
}

}

TEMPLATE_TEST_CASE("PID factory accepts a default and a typical configuration with the first output unchanged",
                   "[pid][config]", float, double)
{
    using Scalar = TestType;
    using Pid = ctrlpp::pid<Scalar, 1>;
    using vector_t = typename Pid::vector_t;
    auto const sp = vector_t::Constant(Scalar{1});
    auto const meas = vector_t::Constant(Scalar{0});
    auto const dt = static_cast<Scalar>(0.01);

    auto idle = ctrlpp::test::constructed(Pid::create(typename Pid::config_type{}));
    auto const idle_output = ctrlpp::test::commanded(idle.compute(sp, meas, dt));
    REQUIRE(same_bits(idle_output, vector_t(vector_t::Zero())));

    typename Pid::config_type cfg{};
    cfg.kp = vector_t::Constant(Scalar{2});
    cfg.ki = vector_t::Constant(Scalar{1});
    cfg.kd = vector_t::Constant(static_cast<Scalar>(0.1));
    cfg.output_min = vector_t::Constant(Scalar{-10});
    cfg.output_max = vector_t::Constant(Scalar{10});
    auto typical = ctrlpp::test::constructed(Pid::create(cfg));
    auto const typical_output = ctrlpp::test::commanded(typical.compute(sp, meas, dt));
    // kp * e + ki * e * dt with e = 1, as the controller rounds it in Scalar.
    if constexpr(std::is_same_v<Scalar, float>)
        REQUIRE(same_bits(typical_output, vector_t(vector_t::Constant(0x1.0147aep+1f))));
    else
        REQUIRE(same_bits(typical_output, vector_t(vector_t::Constant(0x1.0147ae147ae14p+1))));
}

TEMPLATE_TEST_CASE("PID factory names the first fault in check order", "[pid][config]", float, double)
{
    using Scalar = TestType;
    require_refusal_order<filtered_pid<Scalar>>(filtered_baseline<Scalar>(), filtered_faults<Scalar>());
    require_refusal_order<conditional_pid<Scalar>>(conditional_baseline<Scalar>(), conditional_faults<Scalar>());
}

TEMPLATE_TEST_CASE("PID factory refuses each documented condition", "[pid][config]", float, double)
{
    using Scalar = TestType;
    using Pid = parallel_pid<Scalar>;
    using Isa = filtered_pid<Scalar>;
    auto const base = [] {
        auto cfg = with_base_fields(typename Pid::config_type{});
        cfg.template policy<ctrlpp::deriv_filter>().n = {Scalar{10}, Scalar{10}};
        return cfg;
    }();
    auto const isa = filtered_baseline<Scalar>();
    REQUIRE(Pid::create(base).has_value());

    SECTION("gains")
    {
        for(auto const bad : {not_a_number<Scalar>(), infinity<Scalar>(), -infinity<Scalar>()})
        {
            auto kp = base, ki = base, kd = base, ti = isa;
            kp.kp[1] = bad;
            ki.ki[1] = bad;
            kd.kd[1] = bad;
            ti.ki[1] = bad;
            REQUIRE(refusal_of<Pid>(kp) == pid_config_error::non_finite_gain);
            REQUIRE(refusal_of<Pid>(ki) == pid_config_error::non_finite_gain);
            REQUIRE(refusal_of<Pid>(kd) == pid_config_error::non_finite_gain);
            if(bad != infinity<Scalar>())
                REQUIRE(refusal_of<Isa>(ti) == pid_config_error::non_finite_gain);
        }
    }

    SECTION("output limits")
    {
        auto lo_nan = base, hi_nan = base, crossed = base, lo_up = base, hi_down = base;
        lo_nan.output_min[1] = not_a_number<Scalar>();
        hi_nan.output_max[1] = not_a_number<Scalar>();
        crossed.output_min[1] = Scalar{1};
        crossed.output_max[1] = Scalar{0};
        lo_up.output_min[1] = infinity<Scalar>();
        lo_up.output_max[1] = infinity<Scalar>();
        hi_down.output_min[1] = -infinity<Scalar>();
        hi_down.output_max[1] = -infinity<Scalar>();
        for(auto const& cfg : {lo_nan, hi_nan, crossed, lo_up, hi_down})
            REQUIRE(refusal_of<Pid>(cfg) == pid_config_error::invalid_output_limits);

        auto equal = base, unbounded = base;
        equal.output_min[1] = Scalar{1};
        equal.output_max[1] = Scalar{1};
        unbounded.output_min[1] = -infinity<Scalar>();
        unbounded.output_max[1] = infinity<Scalar>();
        REQUIRE(Pid::create(equal).has_value());
        REQUIRE(Pid::create(unbounded).has_value());
    }

    SECTION("setpoint weights")
    {
        for(auto const bad : {not_a_number<Scalar>(), infinity<Scalar>(), -infinity<Scalar>()})
        {
            auto b = base, c = base;
            b.b[1] = bad;
            c.c[1] = bad;
            REQUIRE(refusal_of<Pid>(b) == pid_config_error::non_finite_setpoint_weight);
            REQUIRE(refusal_of<Pid>(c) == pid_config_error::non_finite_setpoint_weight);
        }
    }

    SECTION("filter time constants")
    {
        for(auto const bad : {not_a_number<Scalar>(), infinity<Scalar>(), -infinity<Scalar>(), Scalar{-1}})
        {
            auto sp = base, pv = base;
            sp.template policy<ctrlpp::setpoint_filter>().tf[1] = bad;
            pv.template policy<ctrlpp::pv_filter>().tf[1] = bad;
            REQUIRE(refusal_of<Pid>(sp) == pid_config_error::invalid_filter_time_constant);
            REQUIRE(refusal_of<Pid>(pv) == pid_config_error::invalid_filter_time_constant);
        }
    }

    SECTION("rate limit")
    {
        for(auto const bad : {not_a_number<Scalar>(), -infinity<Scalar>(), Scalar{-1}})
        {
            auto cfg = base;
            cfg.template policy<ctrlpp::rate_limit>().rate_max[1] = bad;
            REQUIRE(refusal_of<Pid>(cfg) == pid_config_error::invalid_rate_limit);
        }
        auto off = base, unbounded = base;
        off.template policy<ctrlpp::rate_limit>().rate_max[1] = Scalar{0};
        unbounded.template policy<ctrlpp::rate_limit>().rate_max[1] = infinity<Scalar>();
        REQUIRE(Pid::create(off).has_value());
        REQUIRE(Pid::create(unbounded).has_value());
    }

    SECTION("anti-windup and derivative filter parameters")
    {
        for(auto const bad : {not_a_number<Scalar>(), infinity<Scalar>(), -infinity<Scalar>()})
        {
            auto kb = base, n = base;
            kb.template policy<ctrlpp::anti_windup<ctrlpp::back_calc>>().kb[1] = bad;
            n.template policy<ctrlpp::deriv_filter>().n[1] = bad;
            REQUIRE(refusal_of<Pid>(kb) == pid_config_error::invalid_back_calculation_gain);
            REQUIRE(refusal_of<Pid>(n) == pid_config_error::invalid_derivative_filter);
        }
        auto threshold = conditional_baseline<Scalar>();
        threshold.template policy<ctrlpp::anti_windup<ctrlpp::conditional_integration>>().error_threshold[1] = not_a_number<Scalar>();
        REQUIRE(refusal_of<conditional_pid<Scalar>>(threshold) == pid_config_error::invalid_error_threshold);
    }

    SECTION("ISA internal gains that overflow")
    {
        auto integral = isa, derivative = isa;
        integral.kp[1] = std::numeric_limits<Scalar>::max();
        integral.ki[1] = static_cast<Scalar>(0.5);
        derivative.kp[1] = std::numeric_limits<Scalar>::max();
        derivative.ki[1] = Scalar{4};
        derivative.kd[1] = Scalar{2};
        REQUIRE(refusal_of<Isa>(integral) == pid_config_error::unrepresentable_gain);
        REQUIRE(refusal_of<Isa>(derivative) == pid_config_error::unrepresentable_gain);
    }

    SECTION("derived filter time constant and default tracking gain that are not representable")
    {
        auto underflow = base, overflow = base, sqrt_kb = base, ratio_kb = base;
        underflow.kd[1] = std::numeric_limits<Scalar>::denorm_min();
        underflow.template policy<ctrlpp::anti_windup<ctrlpp::back_calc>>().kb[1] = Scalar{1};
        overflow.kd[1] = std::numeric_limits<Scalar>::max();
        overflow.template policy<ctrlpp::deriv_filter>().n[1] = static_cast<Scalar>(0.25);
        sqrt_kb.ki[1] = std::numeric_limits<Scalar>::max();
        ratio_kb.ki[1] = std::numeric_limits<Scalar>::max();
        ratio_kb.kp[1] = static_cast<Scalar>(0.25);
        ratio_kb.kd[1] = Scalar{0};
        for(auto const& cfg : {underflow, overflow, sqrt_kb, ratio_kb})
            REQUIRE(refusal_of<Pid>(cfg) == pid_config_error::unrepresentable_gain);
    }
}

TEMPLATE_TEST_CASE("PID set_params refuses without touching the controller", "[pid][config]", float, double)
{
    using Scalar = TestType;
    using Pid = filtered_pid<Scalar>;
    auto const baseline = filtered_baseline<Scalar>();
    auto const sp = Pid::vector_t::Constant(Scalar{3});
    auto const meas = Pid::vector_t::Constant(Scalar{1});
    auto const dt = static_cast<Scalar>(0.01);

    for(auto const& f : filtered_faults<Scalar>())
    {
        auto controller = ctrlpp::test::constructed(Pid::create(baseline));
        auto twin = ctrlpp::test::constructed(Pid::create(baseline));
        drive(controller, 7);
        drive(twin, 7);
        auto bad = baseline;
        f.apply(bad);

        auto const refused = controller.set_params(bad);
        REQUIRE_FALSE(refused.has_value());
        REQUIRE(refused.error() == f.error);
        REQUIRE(same_config(controller.params(), twin.params()));
        REQUIRE(same_policies<Scalar>(controller.params(), twin.params()));
        REQUIRE(same_bits(controller.integral(), twin.integral()));
        auto const next = ctrlpp::test::commanded(controller.compute(sp, meas, dt));
        auto const twin_next = ctrlpp::test::commanded(twin.compute(sp, meas, dt));
        REQUIRE(same_bits(next, twin_next));
    }
}

TEMPLATE_TEST_CASE("PID set_params accepts a valid configuration and keeps the integrator", "[pid][config]", float,
                   double)
{
    using Scalar = TestType;
    using Pid = filtered_pid<Scalar>;
    auto controller = ctrlpp::test::constructed(Pid::create(filtered_baseline<Scalar>()));
    drive(controller, 7);
    auto const integral_before = controller.integral();

    auto retuned = filtered_baseline<Scalar>();
    retuned.kp = Pid::vector_t::Constant(Scalar{3});
    REQUIRE(controller.set_params(retuned).has_value());
    REQUIRE(same_bits(controller.integral(), integral_before));
    REQUIRE(same_bits(controller.params().kp, retuned.kp));
}

TEMPLATE_TEST_CASE("PID set_integral refuses a non-finite value and accepts a finite one", "[pid][config]", float,
                   double)
{
    using Scalar = TestType;
    using Pid = ctrlpp::pid<Scalar, 2>;
    using vector_t = typename Pid::vector_t;
    auto cfg = with_base_fields(typename Pid::config_type{});
    auto controller = ctrlpp::test::constructed(Pid::create(cfg));
    drive(controller, 5);
    auto const integral_before = controller.integral();

    for(auto const bad : {not_a_number<Scalar>(), infinity<Scalar>(), -infinity<Scalar>()})
    {
        for(Eigen::Index channel = 0; channel < 2; ++channel)
        {
            auto seeded = vector_t::Constant(Scalar{1}).eval();
            seeded[channel] = bad;
            auto const refused = controller.set_integral(seeded);
            REQUIRE_FALSE(refused.has_value());
            REQUIRE(refused.error() == pid_config_error::non_finite_integral);
            REQUIRE(same_bits(controller.integral(), integral_before));
        }
    }

    vector_t finite;
    finite << static_cast<Scalar>(0.375), Scalar{-1.5};
    REQUIRE(controller.set_integral(finite).has_value());
    REQUIRE(same_bits(controller.integral(), finite));
}
