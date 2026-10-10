#ifndef HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_SWITCH_POPULATION_H
#define HPP_GUARD_PLANNER_OVERSHOOT_SWEEP_SWITCH_POPULATION_H

#include "common.h"
#include "retarget_windows.h"

#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <cstdint>
#include <algorithm>
#include <type_traits>

namespace planner_overshoot_sweep
{

inline constexpr std::int32_t seeded_pairs = 512;
inline constexpr std::int32_t threshold_steps = 64;
inline constexpr double resolution_multiples[] = {0.25, 0.5, 1.0, 2.0};

enum class switch_leg : std::uint32_t
{
    trapezoidal = 1,
    double_s = 2,
    online_planner_2nd = 3,
    online_planner_3rd = 4,
};

template <typename Scalar>
auto switch_seed(switch_leg leg) -> std::uint32_t
{
    return 100 * static_cast<std::uint32_t>(leg) + (std::is_same_v<Scalar, float> ? 2 : 1);
}

// std::mt19937's output sequence is fixed by the standard and the standard distributions are not, so every value
// is formed from 24 raw bits, which float and double both hold exactly.
template <typename Scalar>
class fraction_source
{
public:
    explicit fraction_source(std::uint32_t seed)
        : engine_{seed}
    {
    }

    auto signed_fraction() -> Scalar { return std::ldexp(static_cast<Scalar>(draw() - (1 << 23)), -23); }

    auto unit_fraction() -> Scalar { return std::ldexp(static_cast<Scalar>(draw()), -24); }

    auto magnitude() -> Scalar { return Scalar{1} + Scalar{63} * unit_fraction(); }

private:
    std::mt19937 engine_;

    auto draw() -> std::int32_t { return static_cast<std::int32_t>(engine_() >> 8); }
};

template <typename Scalar>
struct shape_command
{
    Scalar q0;
    Scalar q1;
    Scalar v0;
    Scalar v1;
    typed_limits<Scalar> limits;
};

struct fraction_pair
{
    double f0;
    double f1;
};

inline constexpr fraction_pair fixed_pairs[] = {{0, 0}, {1, 0}, {0, 1}, {0, 0.5}, {0, -0.5}, {0.5, 0}, {-0.5, 0}};
inline constexpr std::int32_t fixed_pair_count = 7;

template <typename Scalar>
auto next_pair(fraction_source<Scalar>& source, std::int32_t k) -> std::array<Scalar, 2>
{
    if (k < fixed_pair_count) {
        return {static_cast<Scalar>(fixed_pairs[k].f0), static_cast<Scalar>(fixed_pairs[k].f1)};
    }
    auto const f0 = source.signed_fraction();
    return {f0, source.signed_fraction()};
}

template <typename Scalar>
auto acceleration_stretch(typed_limits<Scalar> const& limits) -> Scalar
{
    return limits.v_max / limits.a_max + limits.a_max / limits.j_max;
}

// @cite biagiotti2009 -- Sec. 3.4.1, eq. (3.19)-(3.21): mean speed times the jerk-limited ramp's duration
template <typename Scalar>
auto transition_distance(typed_limits<Scalar> const& limits, Scalar pv_hi, Scalar pv_lo) -> Scalar
{
    auto const rise = pv_hi - pv_lo;
    auto const T = rise * limits.j_max < limits.a_max * limits.a_max ? Scalar{2} * std::sqrt(rise / limits.j_max)
                                                                       : rise / limits.a_max + limits.a_max / limits.j_max;
    return (pv_hi + pv_lo) * T / Scalar{2};
}

template <typename Scalar, typename Visit>
void for_each_trapezoid_command(fraction_source<Scalar>& source, Visit const& visit)
{
    for (auto const& raw : retarget_limits) {
        auto const limits = typed<Scalar>(raw);
        for (auto const sigma : {Scalar{1}, Scalar{-1}}) {
            for (std::int32_t k = 0; k < fixed_pair_count + seeded_pairs; ++k) {
                auto const K = source.magnitude();
                auto const f = next_pair(source, k);
                auto const q1 = sigma * K * limits.v_max * limits.v_max / limits.a_max;
                visit(shape_command<Scalar>{Scalar{0}, q1, sigma * f[0] * limits.v_max, sigma * f[1] * limits.v_max, limits});
            }
        }
    }
}

template <typename Scalar, typename Visit>
void for_each_seeded_double_s_command(fraction_source<Scalar>& source, Visit const& visit)
{
    for (auto const& raw : retarget_limits) {
        auto const limits = typed<Scalar>(raw);
        for (auto const sigma : {Scalar{1}, Scalar{-1}}) {
            for (std::int32_t k = 0; k < fixed_pair_count + seeded_pairs; ++k) {
                auto const K = source.magnitude();
                auto const f = next_pair(source, k);
                auto const pv0 = f[0] * limits.v_max;
                auto const pv1 = f[1] * limits.v_max;
                auto const transition = transition_distance(limits, std::max(pv0, pv1), std::min(pv0, pv1));
                auto const q1 = sigma * K * std::max(transition, limits.v_max * acceleration_stretch(limits));
                visit(shape_command<Scalar>{Scalar{0}, q1, sigma * pv0, sigma * pv1, limits});
            }
        }
    }
}

template <typename Scalar>
auto long_exponents() -> std::array<std::int32_t, 3>
{
    constexpr auto digits = std::numeric_limits<Scalar>::digits;
    return {digits / 4, digits / 2, 3 * digits / 4};
}

template <typename Scalar>
auto resolution_at(Scalar t) -> Scalar
{
    return std::nextafter(t, std::numeric_limits<Scalar>::infinity()) - t;
}

// Durations near 2^m have a time resolution of ulp(2^m), and each jerk phase lasts c of those units.
template <typename Scalar, typename Visit>
void for_each_long_double_s_command(fraction_source<Scalar>& source, Visit const& visit)
{
    for (auto const& raw : retarget_limits) {
        for (auto const sigma : {Scalar{1}, Scalar{-1}}) {
            for (auto const m : long_exponents<Scalar>()) {
                auto const span = std::ldexp(Scalar{1}, m);
                for (auto const c : resolution_multiples) {
                    auto limits = typed<Scalar>(raw);
                    limits.j_max = limits.a_max / (static_cast<Scalar>(c) * resolution_at(span));
                    for (auto const end : {Scalar{-1}, Scalar{1}}) {
                        auto const v0 = sigma * source.signed_fraction() * limits.v_max;
                        visit(shape_command<Scalar>{Scalar{0}, sigma * limits.v_max * span, v0,
                                                    end * sigma * limits.v_max, limits});
                    }
                }
            }
        }
    }
}

// A ramp rising by d is triangular exactly when d j_max < a_max^2, so d is walked across a_max^2 / j_max.
template <typename Scalar, typename Visit>
void for_each_threshold_double_s_command(Visit const& visit)
{
    for (auto const& raw : retarget_limits) {
        auto const base = typed<Scalar>(raw);
        for (auto const sigma : {Scalar{1}, Scalar{-1}}) {
            auto d = base.a_max * base.a_max / base.j_max;
            for (std::int32_t i = 0; i < threshold_steps; ++i) {
                d = std::nextafter(d, Scalar{0});
            }
            for (std::int32_t i = 0; i <= 2 * threshold_steps; ++i, d = std::nextafter(d, std::numeric_limits<Scalar>::infinity())) {
                auto const limits = typed_limits<Scalar>{Scalar{2} * d, base.a_max, base.j_max};
                auto const q1 = sigma * Scalar{1000} * limits.v_max * acceleration_stretch(limits);
                visit(shape_command<Scalar>{Scalar{0}, q1, sigma * d, Scalar{0}, limits});
            }
        }
    }
}

}

#endif
