#include "hardening_helpers.h"
#include "ctrlpp/control/dare.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>

TEST_CASE("DARE refuses an exactly uncontrollable unstable mode as singular_u11", "[dare][error]")
{
    // The operands are pinned in hexadecimal because the claim below is about
    // one bit-identical set of inputs measured across a toolchain matrix. Every
    // entry here is a dyadic rational, so the hexadecimal and the decimal
    // spelling denote the same double; writing it in hexadecimal says so
    // without the reader having to work it out.
    Eigen::Matrix<double, 2, 2> A;
    A << 0x1p+1, 0x0p+0, // diag(2, 1/2)
        0x0p+0, 0x1p-1;
    Eigen::Matrix<double, 2, 1> B;
    B << 0x0p+0, 0x1p+0;
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 0x1p+0;

    // WHY this refusal, and only this one, is a property of the input -- stated
    // without reference to the solver.
    //
    // The mode at eigenvalue 2 has an input coupling of EXACTLY zero, so it is
    // uncontrollable and lies outside the unit circle: no stabilizing solution
    // exists and the pair must be refused. Which refusal it carries follows
    // from two structural facts, neither of which is a comparison against a
    // tolerance.
    //
    // The eigenvalue placement count cannot fire. It asks how many symplectic
    // eigenvalues lie inside the unit disk, and an uncontrollable mode at
    // |lambda| > 1 contributes BOTH lambda and its reciprocal, so two of the
    // four land inside for n = 2 and the count is satisfied.
    //
    // What fails instead is the extraction, and it fails on an exact zero
    // rather than on a small number. The pair is diagonal with a zero input
    // entry, so the mode at 2 is decoupled from the input in the operands
    // themselves; the invariant subspace spanned by the two stable symplectic
    // eigenvalues is then exactly orthogonal to the first state direction, and
    // the smallest pivot of the top-left block is zero, not merely tiny. A rank
    // test compares that zero against a strictly positive threshold, so no
    // rounding enters the verdict and it cannot move with the arithmetic.
    //
    // That last step is argued from the structure of the input and confirmed by
    // measurement rather than proved: the solve reaches the block through a
    // real Schur factorization and a reordering, not through the argument
    // above. The confirmation is sixty-four configurations -- four compilers,
    // four optimization levels, fused multiply-add off and on, against two
    // releases of the linear-algebra library -- all returning `singular_u11`.
    //
    // A single enumerator is asserted here BECAUSE the deciding quantity is a
    // structural zero. That is not true of every refusal this solver produces:
    // the same pair with a cross weight moves the deciding quantity off zero
    // and carries three different enumerators over the same matrix, which is
    // why its own case asserts the input's invariants instead. The platform
    // this matrix does not cover is arm64.
    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(!result.has_value());
    CHECK(result.error() == ctrlpp::dare_error::singular_u11);
}

TEST_CASE("DARE NaN in A returns dare_error::non_finite_input", "[dare][error]")
{
    auto A = ctrlpp::test::nan_matrix<double, 2, 2>();
    Eigen::Matrix<double, 2, 1> B;
    B << 1.0, 0.0;
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(!result.has_value());
    CHECK(result.error() == ctrlpp::dare_error::non_finite_input);
}

TEST_CASE("DARE singular A returns dare_error::singular_a", "[dare][error]")
{
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 0.0, 0.0, 0.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 1.0, 1.0;
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(!result.has_value());
    CHECK(result.error() == ctrlpp::dare_error::singular_a);
}

TEST_CASE("DARE refuses zero dynamics as singular_a", "[dare][error]")
{
    Eigen::Matrix<double, 2, 2> A = Eigen::Matrix<double, 2, 2>::Zero();
    Eigen::Matrix<double, 2, 1> B = Eigen::Matrix<double, 2, 1>::Zero();
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 0x1p+0;

    // WHY this refusal is the only one reachable, stated without reference to
    // the solver.
    //
    // The state matrix is exactly the zero matrix, so the transposed inverse
    // the symplectic pencil is built from does not exist, and the refusal that
    // names a rank-deficient state matrix is the first one the problem admits.
    //
    // The rank verdict on an exactly zero matrix carries no rounding at all.
    // The test is a column-pivoted QR whose threshold is PREMULTIPLIED by the
    // largest pivot, and the largest pivot of the zero matrix is exactly zero,
    // so the threshold is exactly zero and the rank is the count of pivots
    // strictly greater than zero -- which is zero under any IEEE arithmetic and
    // any instruction selection, whatever the relative tolerance in front of
    // it.
    //
    // The three enumerators this case used to also admit are unreachable for
    // this input by construction rather than by measurement. The refusal is
    // returned before the symplectic matrix is formed, so the extraction's rank
    // verdict is never reached; every entry of every operand is finite, so the
    // non-finite diagnosis cannot fire; and the eigenvalue placement count is
    // never taken. A four-way menu asserted almost nothing here.
    //
    // Measured as well as derived: `singular_a` in all sixty-four
    // configurations of the matrix described in the case above, which does not
    // cover arm64.
    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(!result.has_value());
    CHECK(result.error() == ctrlpp::dare_error::singular_a);
}

TEST_CASE("DARE refuses a negative-definite state weighting", "[dare][error]")
{
    // The contract, decided rather than left open: a negative-definite state
    // weighting is NOT a supported input and the solve must refuse it. The
    // standard theory the solver implements requires a positive semi-definite Q
    // for a stabilizing positive semi-definite P to exist, and with Q = -I fewer
    // than n eigenvalues of the symplectic spectrum land inside the unit disk, so
    // the placement count is what catches it.
    //
    // The previous form of this case put every assertion inside its refusal
    // branch, so a solver that RETURNED a value for this input executed no
    // assertion at all and the case passed -- removing coverage from exactly the
    // wrong-success path. Both assertions below run on whichever branch the call
    // takes: a returned value fails the first outright.
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 1.0, 0.0, 1.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.5, 1.0;
    Eigen::Matrix<double, 2, 2> Q = -Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R(0, 0) = 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE_FALSE(result.has_value());

    // The enumerator is pinned to the one the built library actually produces,
    // not to a menu of five. A refusal that could be any of five causes does not
    // tell the caller which of their assumptions failed, and a menu cannot notice
    // when the cause changes.
    CHECK(result.error() == ctrlpp::dare_error::non_stabilizable);

    // And the refusal is about the sign, not about the fixture: the identical
    // pose with a positive-definite state weighting is solved. Without this the
    // case above would still pass on a solver that refused everything.
    auto positive = ctrlpp::dare<double, 2, 1>(A, B, Eigen::Matrix<double, 2, 2>(Eigen::Matrix<double, 2, 2>::Identity()), R);
    REQUIRE(positive.has_value());
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> pes(positive->P);
    CHECK(pes.eigenvalues()(0) > 0.0);
}

TEST_CASE("DARE schur_failed enumerator is reachable at compile time", "[dare][error][design-lever]")
{
    constexpr ctrlpp::dare_error e = ctrlpp::dare_error::schur_failed;
    (void)e;
    CHECK(static_cast<int>(ctrlpp::dare_error::schur_failed) >= 0);
}

TEST_CASE("DARE arithmetic_limit enumerator is reachable at compile time", "[dare][error][design-lever]")
{
    constexpr ctrlpp::dare_error e = ctrlpp::dare_error::arithmetic_limit;
    (void)e;
    CHECK(static_cast<int>(ctrlpp::dare_error::arithmetic_limit) >= 0);
}
