#ifndef HPP_GUARD_CTRLPP_DETAIL_RUNGE_KUTTA_H
#define HPP_GUARD_CTRLPP_DETAIL_RUNGE_KUTTA_H

#include "ctrlpp/types.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>

namespace ctrlpp::detail
{

// ---------------------------------------------------------------------------
// Butcher tableaux for fixed-step Runge-Kutta methods
// ---------------------------------------------------------------------------

struct euler_tableau
{
    static constexpr std::size_t stages = 1;
    static constexpr std::array<double, 1> c{0.0};
    static constexpr std::array<double, 1> b{1.0};
    static constexpr std::array<std::array<double, 1>, 1> a{{{0.0}}};
};

struct rk2_midpoint_tableau
{
    static constexpr std::size_t stages = 2;
    static constexpr std::array<double, 2> c{0.0, 0.5};
    static constexpr std::array<double, 2> b{0.0, 1.0};
    static constexpr std::array<std::array<double, 2>, 2> a{{{0.0, 0.0}, {0.5, 0.0}}};
};

struct rk4_tableau
{
    static constexpr std::size_t stages = 4;
    static constexpr std::array<double, 4> c{0.0, 0.5, 0.5, 1.0};
    static constexpr std::array<double, 4> b{1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0};
    static constexpr std::array<std::array<double, 4>, 4> a{{{0.0, 0.0, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0}, {0.0, 0.5, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}}};
};

// ---------------------------------------------------------------------------
// Traits to extract dimensions from Eigen column vectors
// ---------------------------------------------------------------------------

template <typename T>
struct vector_traits;

template <typename Scalar, int Rows>
struct vector_traits<Eigen::Matrix<Scalar, Rows, 1>>
{
    using scalar_type = Scalar;
    static constexpr std::size_t size = static_cast<std::size_t>(Rows);
};

// ---------------------------------------------------------------------------
// Generic fixed-step Runge-Kutta integrator
// ---------------------------------------------------------------------------

template <typename Tableau, typename F, typename StateVec, typename InputVec>
    requires std::is_same_v<StateVec, std::remove_cvref_t<StateVec>>
auto rk_step(const F& f, const StateVec& x, const InputVec& u, typename vector_traits<StateVec>::scalar_type dt) -> StateVec
{
    using Scalar = typename vector_traits<StateVec>::scalar_type;
    constexpr std::size_t S = Tableau::stages;
    std::array<StateVec, S> k;

    for(std::size_t i = 0; i < S; ++i)
    {
        StateVec x_stage = x;
        for(std::size_t j = 0; j < i; ++j)
        {
            auto aij = static_cast<Scalar>(Tableau::a[i][j]);
            if(aij != Scalar{0})
                x_stage += (dt * aij) * k[j];
        }
        k[i] = f(x_stage, u);
    }

    StateVec x_next = x;
    for(std::size_t i = 0; i < S; ++i)
    {
        auto bi = static_cast<Scalar>(Tableau::b[i]);
        if(bi != Scalar{0})
            x_next += (dt * bi) * k[i];
    }
    return x_next;
}

// ---------------------------------------------------------------------------
// Convenience wrappers
// ---------------------------------------------------------------------------

template <typename F, typename StateVec, typename InputVec>
auto euler_step(const F& f, const StateVec& x, const InputVec& u, typename vector_traits<StateVec>::scalar_type dt) -> StateVec
{
    return rk_step<euler_tableau>(f, x, u, dt);
}

template <typename F, typename StateVec, typename InputVec>
auto rk2_step(const F& f, const StateVec& x, const InputVec& u, typename vector_traits<StateVec>::scalar_type dt) -> StateVec
{
    return rk_step<rk2_midpoint_tableau>(f, x, u, dt);
}

template <typename F, typename StateVec, typename InputVec>
auto rk4_step(const F& f, const StateVec& x, const InputVec& u, typename vector_traits<StateVec>::scalar_type dt) -> StateVec
{
    return rk_step<rk4_tableau>(f, x, u, dt);
}

// ---------------------------------------------------------------------------
// Dormand-Prince RK45 adaptive step-size integrator
// ---------------------------------------------------------------------------

template <typename Scalar>
struct rk45_config
{
    Scalar atol = Scalar{1e-6};
    Scalar rtol = Scalar{1e-6};
    Scalar dt_min = Scalar{1e-12};
    Scalar dt_max = Scalar{1.0};
};

template <typename Scalar, std::size_t NX>
struct rk45_result
{
    Vector<Scalar, NX> x_next;
    Scalar dt_actual;
    Scalar dt_next;
};

template <typename F, typename StateVec, typename InputVec>
auto rk45_step(const F& f, const StateVec& x, const InputVec& u, typename vector_traits<StateVec>::scalar_type dt, const rk45_config<typename vector_traits<StateVec>::scalar_type>& cfg = {})
    -> rk45_result<typename vector_traits<StateVec>::scalar_type, vector_traits<StateVec>::size>
{
    using Scalar = typename vector_traits<StateVec>::scalar_type;
    constexpr std::size_t NX = vector_traits<StateVec>::size;

    // Dormand-Prince coefficients
    static constexpr double a21 = 1.0 / 5.0;
    static constexpr double a31 = 3.0 / 40.0, a32 = 9.0 / 40.0;
    static constexpr double a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0;
    static constexpr double a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0;
    static constexpr double a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0;

    // 5th order weights
    static constexpr double b1 = 35.0 / 384.0, b3 = 500.0 / 1113.0, b4 = 125.0 / 192.0, b5 = -2187.0 / 6784.0, b6 = 11.0 / 84.0;

    // Error weights: e_i = b_i - b*_i (difference between 5th and 4th order)
    static constexpr double e1 = 71.0 / 57600.0, e3 = -71.0 / 16695.0, e4 = 71.0 / 1920.0, e5 = -17253.0 / 339200.0, e6 = 22.0 / 525.0, e7 = -1.0 / 40.0;

    Scalar h = std::clamp(dt, cfg.dt_min, cfg.dt_max);

    for(;;)
    {
        StateVec k1 = f(x, u);
        StateVec k2 = f(StateVec{(x + h * Scalar(a21) * k1).eval()}, u);
        StateVec k3 = f(StateVec{(x + h * (Scalar(a31) * k1 + Scalar(a32) * k2)).eval()}, u);
        StateVec k4 = f(StateVec{(x + h * (Scalar(a41) * k1 + Scalar(a42) * k2 + Scalar(a43) * k3)).eval()}, u);
        StateVec k5 = f(StateVec{(x + h * (Scalar(a51) * k1 + Scalar(a52) * k2 + Scalar(a53) * k3 + Scalar(a54) * k4)).eval()}, u);
        StateVec k6 = f(StateVec{(x + h * (Scalar(a61) * k1 + Scalar(a62) * k2 + Scalar(a63) * k3 + Scalar(a64) * k4 + Scalar(a65) * k5)).eval()}, u);

        // 5th order solution
        StateVec x5 = (x + h * (Scalar(b1) * k1 + Scalar(b3) * k3 + Scalar(b4) * k4 + Scalar(b5) * k5 + Scalar(b6) * k6)).eval();

        StateVec k7 = f(x5, u);

        // Error estimate
        StateVec err_vec = (h * (Scalar(e1) * k1 + Scalar(e3) * k3 + Scalar(e4) * k4 + Scalar(e5) * k5 + Scalar(e6) * k6 + Scalar(e7) * k7)).eval();

        // Compute scaled error norm (mixed absolute/relative tolerance)
        Scalar err{0};
        for(std::size_t i = 0; i < NX; ++i)
        {
            Scalar sc = cfg.atol + cfg.rtol * std::abs(x(static_cast<int>(i)));
            Scalar ei = err_vec(static_cast<int>(i)) / sc;
            err += ei * ei;
        }
        err = std::sqrt(err / Scalar(NX));

        if(err <= Scalar{1} || h <= cfg.dt_min)
        {
            // Accept step
            Scalar h_new = (err > Scalar{0}) ? h * std::clamp(Scalar{0.84} * std::pow(Scalar{1} / err, Scalar{0.2}), Scalar{0.1}, Scalar{4}) : h * Scalar{4};
            h_new = std::clamp(h_new, cfg.dt_min, cfg.dt_max);

            return {x5, h, h_new};
        }

        // Reject step: shrink h and retry
        Scalar h_new = h * std::max(Scalar{0.84} * std::pow(Scalar{1} / err, Scalar{0.2}), Scalar{0.1});
        h = std::max(h_new, cfg.dt_min);
    }
}

} // namespace ctrlpp::detail

#endif
