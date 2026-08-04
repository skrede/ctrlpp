/// Which state dimensions the discrete solver can actually answer in `float`.
///
/// WHY THIS EXISTS. An on-target capture found that at `float` this damped-chain
/// family is refused at several dimensions while `double` answers every one of
/// them. That is a property of the LIBRARY and of the scalar type, not of the
/// board: it reproduces exactly on the host. It was published in the real-time
/// matrix as current behavior with nothing pinning it, which is the defect this
/// round spent a finding on -- so this file is that finding applied to its own
/// output.
///
/// THE ACCEPTANCE QUANTITY is the estimated relative forward error against the
/// half-significand margin `sqrt(eps)`, which for `float` is 3.4527e-04. The
/// gate asks that half of float's 24-bit significand survive the solve.
///
/// THE STATE DIMENSION IS NOT WHAT DECIDES THIS. Varying the input count
/// independently shows `NX = 8` with four inputs ACCEPTED (0.725 of the margin)
/// while `NX = 6` with three is refused (1.274). What drives the error is
/// `group = NX / NU`, the number of states each input must reach through:
/// `group = 2` lands near `3e-04`, `group = 3` near `9e-04`, `group = 4` to `5`
/// near `1.5e-03`. The margin is `3.4527e-04`, so the whole `group = 2` family
/// straddles it and the accept-or-refuse outcome there is rounding noise.
///
/// The estimator is not the problem: measured against a `double` reference the
/// estimate matches the true relative error to three significant figures at
/// every pose, so the gate refuses answers that really are that wrong.
///
/// This file therefore asserts only the `group >= 3` band, where the error
/// exceeds the margin by a factor of two and a half or more and the verdict is
/// not in doubt. See the real-time matrix for the full sweep.

#include "ctrlpp/control/dare.h"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <algorithm>

namespace
{

/// The discrete damped chain the benchmarks and the on-target probe both sweep:
/// a forward-Euler step of a chain with -0.5 on the diagonal and 1.0 on the
/// superdiagonal. `Q = I` and `R = 0.1 I` put the weight scale at exactly one,
/// so the entry point runs its acceptance check ONCE and this file measures that
/// check rather than the equilibrated path's two.
template <typename Scalar, std::size_t NX, std::size_t NU>
struct damped_chain
{
    static constexpr int n  = int(NX);
    static constexpr int nu = int(NU);

    Eigen::Matrix<Scalar, n, n>   A = Eigen::Matrix<Scalar, n, n>::Identity();
    Eigen::Matrix<Scalar, n, nu>  B = Eigen::Matrix<Scalar, n, nu>::Zero();
    Eigen::Matrix<Scalar, n, n>   Q = Eigen::Matrix<Scalar, n, n>::Identity();
    Eigen::Matrix<Scalar, nu, nu> R = Scalar(0.1) * Eigen::Matrix<Scalar, nu, nu>::Identity();

    damped_chain()
    {
        const Scalar dt = Scalar(0.01);
        for(std::size_t i = 0; i < NX; ++i)
            A(int(i), int(i)) += dt * Scalar(-0.5);
        for(std::size_t i = 0; i + 1 < NX; ++i)
            A(int(i), int(i + 1)) = dt;
        const std::size_t group = NX / NU;
        for(std::size_t j = 0; j < NU; ++j)
            B(int((j + 1) * group - 1), int(j)) = dt;
    }
};

/// The acceptance quantity itself, as a fraction of the margin it is compared
/// against. Formed from the solver's OWN answer, the way the solver forms it, so
/// the number this file asserts on is the number the acceptance rule sees.
template <typename Scalar, std::size_t NX, std::size_t NU>
auto forward_error_over_margin() -> Scalar
{
    constexpr int n  = int(NX);
    constexpr int nu = int(NU);
    const damped_chain<Scalar, NX, NU> pose;

    auto operands = ctrlpp::detail::factor_dare_symplectic_operands<Scalar, NX, NU>(
        pose.A, pose.B, pose.R);
    REQUIRE(operands.has_value());
    auto Z = ctrlpp::detail::build_dare_symplectic<Scalar, NX>(pose.A, pose.Q, *operands);
    REQUIRE(Z.has_value());
    auto solved = ctrlpp::detail::dare_solve_from_symplectic<Scalar, NX, NU>(*Z);
    // The SOLVE succeeds at every dimension here. What differs is whether the
    // answer it produced survives the acceptance check, and that distinction is
    // the whole point of the file.
    REQUIRE(solved.has_value());

    const Eigen::Matrix<Scalar, n, n> P = solved->P;
    Eigen::Matrix<Scalar, nu, n>      K;
    REQUIRE(ctrlpp::detail::compute_dare_gain<Scalar, NX, NU>(pose.A, pose.B, pose.R, P, K)
            == ctrlpp::detail::dare_verification::verified);

    const Eigen::Matrix<Scalar, n, n> residual =
        (pose.A.transpose() * P * pose.A - P
         - pose.A.transpose() * P * pose.B * K + pose.Q).eval();
    const Eigen::Matrix<Scalar, n, n> closed_loop = (pose.A - pose.B * K).eval();

    const auto estimate =
        ctrlpp::detail::estimate_riccati_forward_error<Scalar, n>(closed_loop, residual, P);
    REQUIRE(estimate.resolved);
    return estimate.value / ctrlpp::detail::half_significand_margin<Scalar>();
}

}

TEST_CASE("float refuses the damped chain by chain-length-per-input, and double does not")
{
    // THE REFUSED BAND: every pose here has `group >= 3`. Each exceeds the margin
    // by a factor of at least two and a half, so which side of it they land on
    // does not move with instruction selection. The bound is deliberately 2.0
    // rather than the measured 2.590, so a change that halves the estimate still
    // fails this case instead of sliding under a fitted threshold.
    CHECK(forward_error_over_margin<float, 3, 1>() > 2.0F);
    CHECK(forward_error_over_margin<float, 5, 1>() > 2.0F);
    CHECK(forward_error_over_margin<float, 6, 2>() > 2.0F);
    CHECK(forward_error_over_margin<float, 8, 2>() > 2.0F);

    // and the entry point declines them, through the shared enumerator
    CHECK_FALSE(ctrlpp::dare<float, 3, 1>(damped_chain<float, 3, 1>{}.A, damped_chain<float, 3, 1>{}.B,
                                          damped_chain<float, 3, 1>{}.Q, damped_chain<float, 3, 1>{}.R)
                    .has_value());
    CHECK_FALSE(ctrlpp::dare<float, 6, 2>(damped_chain<float, 6, 2>{}.A, damped_chain<float, 6, 2>{}.B,
                                          damped_chain<float, 6, 2>{}.Q, damped_chain<float, 6, 2>{}.R)
                    .has_value());

    // THE COMFORTABLE END. `NX = 2` clears the margin by a factor of nearly six
    // and is the only pose in the sweep that does; every other `group = 2` pose
    // straddles it.
    CHECK(forward_error_over_margin<float, 2, 1>() < 0.5F);
    {
        const damped_chain<float, 2, 1> pose;
        CHECK(ctrlpp::dare<float, 2, 1>(pose.A, pose.B, pose.Q, pose.R).has_value());
    }

    // `NX = 4` IS NOT ASSERTED EITHER WAY, and neither is any other `group = 2`
    // pose. It measured 0.966 of the margin, and
    // three and a half percent is inside what a different instruction selection
    // moves. Recording the bound that IS safe rather than the outcome that is
    // not: whatever side it lands on, it is far below the refused band, so the
    // ordering of the two bands is what this pins.
    CHECK(forward_error_over_margin<float, 4, 2>() < 2.0F);

    // THE MECHANISM, pinned as an ORDERING rather than as any single verdict.
    //
    // What separates the bands is `group = NX / NU`, the number of states each
    // input must reach through -- NOT the state dimension. `NX = 8` with four
    // inputs is a `group = 2` pose and sits BELOW `NX = 6` with two inputs, which
    // is `group = 3`, even though it is the larger problem. Asserting the
    // ordering pins that mechanism while asserting nothing about which side of
    // the margin the `group = 2` family lands on, because it straddles it.
    const float g2_worst = std::max({forward_error_over_margin<float, 2, 1>(),
                                     forward_error_over_margin<float, 4, 2>(),
                                     forward_error_over_margin<float, 8, 4>(),
                                     forward_error_over_margin<float, 12, 6>()});
    const float g3plus_best = std::min({forward_error_over_margin<float, 3, 1>(),
                                        forward_error_over_margin<float, 6, 2>(),
                                        forward_error_over_margin<float, 8, 2>(),
                                        forward_error_over_margin<float, 5, 1>()});
    CHECK(g2_worst < g3plus_best);

    // DOUBLE ANSWERS EVERY ONE OF THEM. The refusals above are a property of the
    // scalar type, not of the family or of the solver's construction.
    CHECK(forward_error_over_margin<double, 3, 1>() < 0.5);
    CHECK(forward_error_over_margin<double, 5, 1>() < 0.5);
    CHECK(forward_error_over_margin<double, 6, 2>() < 0.5);
    CHECK(forward_error_over_margin<double, 8, 2>() < 0.5);
    {
        const damped_chain<double, 8, 2> pose;
        CHECK(ctrlpp::dare<double, 8, 2>(pose.A, pose.B, pose.Q, pose.R).has_value());
    }
}
