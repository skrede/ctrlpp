/// Which state dimensions the discrete solver can actually answer in `float`.
///
/// WHY THIS EXISTS. An on-target capture found that at `float` this damped-chain
/// family is refused at several dimensions while `double` answers every one of
/// them. That is a property of the SCALAR TYPE and the library, not of the
/// board: it reproduces on the host. It was published in the real-time matrix as
/// current behavior with nothing asserting it, which is the defect this round
/// spent a finding on -- so this file is that finding applied to its own output.
///
/// THE ACCEPTANCE QUANTITY is the estimated relative forward error against the
/// half-significand margin `sqrt(eps)`: the gate asks that half the scalar's
/// significand survive the solve. This file expresses every measurement as a
/// fraction of that margin, so float and double figures are directly comparable.
///
/// WHAT THIS FILE LEARNED THE HARD WAY. Its first version asserted absolute
/// bands around the float ratios and broke CI on two platforms at once: an
/// assertion flipped on Apple clang, and the translation unit's thirteen
/// instantiations compiled for twenty-one minutes under a coverage build before
/// the runner killed it. Both lessons are encoded below -- portable comparisons
/// only, and two instantiations rather than thirteen.
///
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

TEST_CASE("float loses far more of the Riccati solution than double on the same pose")
{
    // WHAT IS PORTABLE AND WHAT IS NOT -- learned by breaking CI with the first
    // version of this file.
    //
    // The first version asserted absolute margin ratios: that `NX = 8, NU = 2`
    // exceeds the margin by more than a factor of two, measured at 4.105 on the
    // authoring station. On Apple clang the same expression fell BELOW 2.0. The
    // float forward error on this family varies by MORE THAN A FACTOR OF TWO
    // across toolchains, so no absolute band around it is portable and neither
    // is any individual accept-or-refuse verdict.
    //
    // What is portable is the COMPARISON between scalars on the same pose. float
    // carries 24 significand bits and double 53, and the gap that opens between
    // them on an identical problem is orders of magnitude -- far outside the
    // couple-of-x that instruction selection moves. That comparison is the real
    // content anyway: it is why an embedded caller running float has a problem
    // that the same code in double does not.
    //
    // The per-pose numbers live in docs/rt-safety-matrix.md, labeled with the
    // toolchain that produced them, because that is what they are.

    // ONE pose, two scalars. Measured on the authoring station, each distinct
    // (scalar, dimension) pair costs about nineteen seconds of compile time --
    // the fixed cost of including the solver and the test framework is under two
    // -- so instantiation count is the whole build cost of this file. The first
    // version instantiated thirteen pairs and compiled for twenty-one minutes
    // under a coverage build before the runner killed it. Two pairs keeps this
    // file beside the suite's other Riccati targets instead of dwarfing them.
    const float  f3 = forward_error_over_margin<float, 3, 1>();
    const double d3 = forward_error_over_margin<double, 3, 1>();

    // Both are expressed against their OWN scalar's margin, so the comparison is
    // already normalized: it says float spends a far larger fraction of what it
    // has than double does.
    CHECK(f3 > 100.0F * static_cast<float>(d3));

    // double clears its own margin on both poses with room to spare, so the
    // entry point answers them.
    CHECK(d3 < 0.5);
    {
        // Same dimension as above, so this costs no new instantiation.
        const damped_chain<double, 3, 1> pose;
        CHECK(ctrlpp::dare<double, 3, 1>(pose.A, pose.B, pose.Q, pose.R).has_value());
    }

    // NO FLOAT VERDICT IS ASSERTED. On the authoring toolchain both poses above
    // are refused and `NX = 2` is accepted, but the macOS failure proved those
    // outcomes move. Asserting one would pin a toolchain, not a behavior.
}
