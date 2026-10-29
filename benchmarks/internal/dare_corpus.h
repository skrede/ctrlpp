#ifndef HPP_GUARD_BENCHMARKS_INTERNAL_DARE_CORPUS_H
#define HPP_GUARD_BENCHMARKS_INTERNAL_DARE_CORPUS_H

#include <Eigen/Dense>

#include <tuple>
#include <cstddef>
#include <algorithm>

namespace ctrlpp::bench
{

constexpr char const* dare_residual_metric =
    "relative residual of this arm's own discrete Riccati solution";

/// Forward-Euler step of the continuous damped chain the CARE bakeoff uses:
/// -0.5 on the diagonal and 1.0 on the superdiagonal, one input per group of
/// states. The step keeps every eigenvalue of A at 1 - 0.5 dt, inside the unit
/// disk, so the discrete equation is well posed at every size.
template <std::size_t NX, std::size_t NU>
auto build_discrete_damped_chain()
{
    constexpr int n  = int(NX);
    constexpr int nu = int(NU);
    constexpr double dt = 0.01;

    Eigen::Matrix<double, n, n> A = Eigen::Matrix<double, n, n>::Identity();
    for (std::size_t i = 0; i < NX; ++i)
        A(int(i), int(i)) += dt * -0.5;
    for (std::size_t i = 0; i + 1 < NX; ++i)
        A(int(i), int(i + 1)) = dt * 1.0;

    Eigen::Matrix<double, n, nu> B = Eigen::Matrix<double, n, nu>::Zero();
    const std::size_t group = NX / NU;
    for (std::size_t j = 0; j < NU; ++j)
    {
        const std::size_t last_row = std::min((j + 1) * group, NX) - 1;
        B(int(last_row), int(j)) = dt;
    }

    Eigen::Matrix<double, n, n>   Q = Eigen::Matrix<double, n, n>::Identity();
    Eigen::Matrix<double, nu, nu> R = 0.1 * Eigen::Matrix<double, nu, nu>::Identity();

    return std::tuple{A, B, Q, R};
}

/// The same chain at a weight scale of `s` instead of one.
///
/// `dare_weight_scale` returns the smallest power of two at or above
/// max(||Q||_max, ||R||_max). For `s` a power of two the largest entry is `s`
/// itself, so it returns exactly `s`, and the entry point's gate on that value
/// being different from one opens.
///
/// BOTH WEIGHTS ARE MULTIPLIED AND NOT ONLY Q, deliberately. Scaling Q alone
/// moves the weight ratio with the scale, so the equilibrated problem the solve
/// actually sees would be a different problem at every rung and the difference
/// against the scale-one row would confound the equilibrated path's cost with a
/// change in the work underneath it. Multiplying and dividing by a power of two
/// are exact in the normal range, so Q/s is the identity and R/s is fl(0.1)
/// bitwise: the Schur solve underneath is the scale-one solve at every rung and
/// the difference is the equilibrated path plus the two rescaling passes.
template <std::size_t NX, std::size_t NU>
auto build_scaled_damped_chain(double s)
{
    auto [A, B, Q, R] = build_discrete_damped_chain<NX, NU>();
    Q *= s;
    R *= s;
    return std::tuple{A, B, Q, R};
}

struct weight_scale_rung
{
    double      scale;
    char const* label;
};

/// A geometric ladder of powers of two in the weight scale: the power of two
/// nearest in exponent to each of 1e-6, 1e-4, 1e-2, 1e2, 1e4 and 1e6, plus one
/// rung immediately beside one.
///
/// The overhead is not scale-independent by construction: the equilibrated
/// operands are the scale-one corpus at every rung, but the rescaled answer,
/// its definiteness factorization, the caller-scale acceptance check and the
/// gain agreement run on operands whose exponents move with the scale. A ladder
/// answers whether that cost depends on the scale; one value would answer
/// whether it depends on that value.
///
/// The rung at 2 separates two different statements -- "the gate is open" and
/// "the pose is far from equilibrated" -- because everything below the gate
/// runs there exactly as it runs at 2^20. It is the smallest power of two
/// above one.
constexpr weight_scale_rung scale_rungs[] = {
    {0x1p-20, "2^-20"},
    {0x1p-13, "2^-13"},
    {0x1p-7,  "2^-7" },
    {0x1p+1,  "2"    },
    {0x1p+7,  "2^7"  },
    {0x1p+13, "2^13" },
    {0x1p+20, "2^20" },
};

/// The fourth quantity, one row per rung: the same entry point on weights whose
/// scale is not one, so the call carries the equilibrated path on top of
/// everything the first row already measures. Subtract the scaled row from the
/// scale-one row at the same size to obtain that path's cost.


/// The size sweep, written once. Both the timing set and the corpus check
/// traverse it, so a confirmation cannot silently cover a different set of
/// sizes than the measurement it licenses.
///
/// NX = 15 IS THE CEILING, AND IT IS A COMPILE-TIME ONE. The forward-error
/// estimate holds an M-by-M operator with M = NX(NX+1)/2 as a fixed-size Eigen
/// object, so Eigen's 128 KiB stack-allocation limit admits M <= 128, that is
/// NX <= 15. At NX = 16 the operator is 136 x 136 = 147,968 bytes and the static
/// assertion fires: the acceptance check cannot be instantiated there at all, so
/// there is no larger size for this sweep to reach.
template <typename Fn>
void for_each_size(Fn&& fn)
{
    fn.template operator()<2, 1>();
    fn.template operator()<4, 2>();
    fn.template operator()<6, 2>();
    fn.template operator()<8, 2>();
    fn.template operator()<12, 3>();
    fn.template operator()<15, 3>();
}

}

#endif
