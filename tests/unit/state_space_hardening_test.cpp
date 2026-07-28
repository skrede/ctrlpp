// What the oracles in this file decide.
//
// The discrete state-space type is a plain container for the four matrices of a
// linear recursion; it performs no solve and has no failure channel. What these
// cases can therefore decide is exactly what the recursion's own algebra fixes,
// and the file asserts that and nothing weaker:
//
//  * Exact equality wherever the answer is representable and the arithmetic
//    necessarily reaches it -- an all-zero system's outputs, an identity
//    system's pass-through, and the bitwise-unchanged coordinate of a decoupled
//    recursion. A tolerance on an exact quantity is a weaker statement than the
//    code supports.
//  * A one-addition budget where a decimal literal and the sum that produces it
//    are two roundings of the same real number and so may differ by one ulp.
//  * The specific outcome, not a disjunction, where a non-finite matrix entry
//    enters: NaN in, NaN out, named.
//
// What they deliberately do not decide. Nothing here asserts anything about
// conditioning: the near-singular system below is exercised because its
// recursion is diagonal and therefore analytically closed, not because a
// condition number is being measured. And nothing here validates the matrices
// at construction, because the type does not: it stores what it is given.

#include "hardening_helpers.h"

#include "ctrlpp/model/state_space.h"
#include "ctrlpp/types.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

constexpr double ss_eps = std::numeric_limits<double>::epsilon();

// One step of the identity recursion is a single addition per component: the
// products against the identity are exact and the products against the zero
// feedthrough contribute exact zeros, so the computed sum can differ from the
// decimal literal naming the same real number by at most one rounding at the
// scale of the larger operand.
constexpr int identity_step_ops = 1;

// Rounded operations reaching the fixed point of the decoupled recursion. The
// map is one multiply and one addition per step, and the iteration is a
// contraction that reaches its floating-point fixed point in two steps, so two
// steps of two operations bound the iterate; the closed form it is compared
// against costs one subtraction and one division. Six in all.
constexpr int fixed_point_ops = 6;

}

// ── State space hardening ──────────────────────────────────────────────────────

TEST_CASE("State space with all-zero matrices", "[state_space][hardening][negative]")
{
    auto A = ctrlpp::test::zero_matrix<double, 2, 2>();
    auto B = ctrlpp::test::zero_matrix<double, 2, 1>();
    auto C = ctrlpp::test::zero_matrix<double, 1, 2>();
    auto D = ctrlpp::test::zero_matrix<double, 1, 1>();

    ctrlpp::discrete_state_space<double, 2, 1, 1> sys{A, B, C, D};

    // Propagate: x(k+1) = A*x(k) + B*u(k) = 0
    ctrlpp::Vector<double, 2> x = ctrlpp::Vector<double, 2>::Ones();
    ctrlpp::Vector<double, 1> u;
    u << 1.0;

    // Every product is an exact zero and every sum of exact zeros is an exact
    // zero, so no rounding can enter and the answer is the zero vector itself.
    // The tolerance this replaces admitted a system that leaked its state.
    auto x_next = (sys.A * x + sys.B * u).eval();
    REQUIRE(x_next == ctrlpp::Vector<double, 2>::Zero());

    auto y = (sys.C * x + sys.D * u).eval();
    REQUIRE(y(0) == 0.0);
}

TEST_CASE("State space with NaN in system matrices", "[state_space][hardening][negative]")
{
    auto A = ctrlpp::test::nan_matrix<double, 2, 2>();
    auto B = ctrlpp::Matrix<double, 2, 1>::Identity();
    auto C = ctrlpp::Matrix<double, 1, 2>::Identity();
    auto D = ctrlpp::test::zero_matrix<double, 1, 1>();

    ctrlpp::discrete_state_space<double, 2, 1, 1> sys{A, B, C, D};

    ctrlpp::Vector<double, 2> x = ctrlpp::Vector<double, 2>::Ones();
    ctrlpp::Vector<double, 1> u;
    u << 1.0;

    auto x_next = (sys.A * x + sys.B * u).eval();
    CHECK(std::isnan(x_next(0)));
    CHECK(std::isnan(x_next(1)));
}

TEST_CASE("Identity system propagation matches analytical", "[state_space][hardening][precision]")
{
    // A=I, B=I, C=I, D=0 -> x(k+1) = x(k) + u(k), y(k) = x(k)
    auto A = ctrlpp::Matrix<double, 2, 2>::Identity();
    auto B = ctrlpp::Matrix<double, 2, 2>::Identity();
    auto C = ctrlpp::Matrix<double, 2, 2>::Identity();
    auto D = ctrlpp::test::zero_matrix<double, 2, 2>();

    ctrlpp::discrete_state_space<double, 2, 2, 2> sys{A, B, C, D};

    ctrlpp::Vector<double, 2> x;
    x << 1.0, 2.0;

    ctrlpp::Vector<double, 2> u;
    u << 0.5, -0.3;

    // Step 1: x_next = x + u = [1.5, 1.7]
    auto x_next = (sys.A * x + sys.B * u).eval();
    // 1.0 + 0.5 is exact in binary, so the first component admits no tolerance.
    REQUIRE(x_next(0) == 1.5);
    // 1.7 is not representable, and neither is 0.3. The computed sum and the
    // decimal literal are two roundings of the same real number, so they agree
    // to within one addition at the scale of the larger operand -- which is a
    // statement about the arithmetic, not a guess about the answer.
    REQUIRE_THAT(x_next(1), WithinAbs(1.7, identity_step_ops * ss_eps * 2.0));

    // y = C x with C the identity and D zero, so each output is its own state
    // component bitwise: the products are exact and the added feedthrough is an
    // exact zero.
    auto y = (sys.C * x + sys.D * u).eval();
    REQUIRE(y(0) == 1.0);
    REQUIRE(y(1) == 2.0);

    // Step 2: x_next2 = x_next + u = [2.0, 1.4]
    auto x_next2 = (sys.A * x_next + sys.B * u).eval();
    REQUIRE(x_next2(0) == 2.0);
    // Two accumulated additions now separate the iterate from the literal.
    REQUIRE_THAT(x_next2(1), WithinAbs(1.4, 2 * identity_step_ops * ss_eps * 2.0));
}

TEST_CASE("Near-singular system converges to its closed-form fixed point",
          "[state_space][hardening][robustness]")
{
    constexpr double a = 1e-10;
    auto sys = ctrlpp::test::near_singular_system<double, 2, 1, 1>(a);

    ctrlpp::Vector<double, 2> x;
    x << 1.0, 1.0;

    ctrlpp::Vector<double, 1> u;
    u << 1.0;

    // The builder gives A = diag(a, 1) and B = [1; 0], so the recursion is
    // DIAGONAL and closed in both coordinates:
    //   x0(k+1) = a*x0(k) + 1   -- a contraction with ratio a, whose fixed point
    //                             is 1/(1 - a) and which reaches it in two steps
    //                             because a^2 is already below one ulp there;
    //   x1(k+1) = x1(k)         -- driven by nothing, so bitwise unchanged.
    for(int i = 0; i < 100; ++i)
    {
        x = (sys.A * x + sys.B * u).eval();
    }

    // Bitwise, not within a tolerance: the second row of A is the identity row
    // and the second row of B is exactly zero, so every step multiplies by one
    // and adds an exact zero. A recursion that leaked the first coordinate into
    // the second fails this and would pass any tolerance loose enough to cover
    // the first coordinate's rounding.
    REQUIRE(x(1) == 1.0);

    // The first state does NOT decay: it is driven by B*u = 1 every step and
    // converges UPWARD to 1/(1 - a) = 1.0000000001. The bound this replaces
    // ("bounded, not growing", |x0| < 2) passed for that reason and for the
    // opposite one equally, and the comment beside it was false.
    REQUIRE_THAT(x(0), WithinRel(1.0 / (1.0 - a), fixed_point_ops * ss_eps));
}
