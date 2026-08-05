/// @file
/// @brief Probe: what decides the discrete Riccati solver's error enumerator
/// below its own resolution boundary -- the input, or instruction selection?
///
/// The continuous solver's enumerator was measured to be decided by instruction
/// selection there: one bit-identical input yields three different enumerators
/// across three toolchains. The discrete twin plausibly has the same property
/// and had never been swept. This program is the instrument that settles it; it
/// is driven across compilers, optimization levels and contraction settings by
/// `dare_error_toolchain_sweep.sh` beside it.
///
/// The sweep axis is the one-parameter family the enumerator's own
/// documentation already records as measured:
///
///     A = diag(2, 1/2),  B = [delta; 1],  Q = I,  R = 1
///
/// which is controllable and detectable for every delta != 0, and whose
/// accept/refuse transition that documentation places between 1e-7 and 1e-8.
/// Starting there means the measurement names the deciding quantity rather than
/// being fitted to a quantity chosen in advance.
///
/// Every entry of the family that does not move with the sweep is written as a
/// hexadecimal float literal, and the swept parameter is either a hexadecimal
/// literal or an exact power of two, so a point is bit-identical under any
/// conforming compiler. A decimal literal does not survive the round trip a
/// toolchain table depends on.
///
/// **This program compares against no tolerance of its own and reports no
/// verdict.** It prints what the solver returned at each input coordinate and
/// what the candidate deciding quantities were there. A threshold inside the
/// instrument would be an underived constant of exactly the kind this
/// measurement exists to remove.
///
/// Three quantities are printed at every point, and which of them tracks the
/// transition is left to the reader of the output:
///
///   * visibility  -- how visible the mode nearest the unit circle is through
///                    the state weighting, computed from the weighting's FACTOR.
///                    Computing it from the assembled weighting cancels the
///                    information away and reports exactly zero.
///   * pivot_ratio -- the rank-revealing pivot ratio of the leading block of an
///                    orthonormal basis of the stable invariant subspace. This
///                    is the quantity the solver's extraction compares against
///                    its rank threshold, reconstructed here by a route the
///                    solver does not take.
///   * controllability -- the coupling of the mode outside the unit circle into
///                    the input, |w* B| for the left eigenvector w of that mode.
///                    Printed because the claim "delta IS that coupling" is a
///                    property of A that should be measured rather than assumed.
///
/// Building: standalone, no test framework and no build system. Eigen enters as
/// a system include, the way the library's own build treats it, so the only
/// diagnostics a build reports are this file's own.
///
///     g++ -std=c++20 -O2 -fno-exceptions -fno-rtti
///         -I lib/ctrlpp/include -isystem /usr/include/eigen3
///         tools/dare_error_toolchain_sweep.cpp -o /tmp/dare_error_probe
///
/// Usage: no argument prints the full report; `--anchor` prints the single line
/// the driver assembles its table from.

#include "ctrlpp/control/dare.h"

#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <bit>
#include <cmath>
#include <cstdio>
#include <limits>
#include <complex>
#include <cstdint>
#include <cstring>

namespace {

constexpr int nx = 2;
constexpr int nu = 1;

using state_matrix      = Eigen::Matrix<double, nx, nx>;
using input_matrix      = Eigen::Matrix<double, nx, nu>;
using weight_matrix     = Eigen::Matrix<double, nx, nx>;
using input_weight      = Eigen::Matrix<double, nu, nu>;
using symplectic_matrix = Eigen::Matrix<double, 2 * nx, 2 * nx>;
using subspace_basis    = Eigen::Matrix<double, 2 * nx, nx>;
using column_vector     = Eigen::Matrix<double, nx, 1>;

/// The plain object Eigen deduces for the transpose of a square column-major
/// matrix. Naming it lets the factorizations below take the transposed operand
/// directly, without an intermediate copy for the factorization to bind to.
using transposed_state_matrix = Eigen::Matrix<double, nx, nx, Eigen::RowMajor>;

/// The unstable and stable modes of the family, pinned exactly.
auto family_state() -> state_matrix
{
    state_matrix A;
    A << 0x1p+1, 0x0p+0,
         0x0p+0, 0x1p-1;
    return A;
}

/// The input map. Only the swept entry moves; the other is pinned exactly.
auto family_input(double delta) -> input_matrix
{
    input_matrix B;
    B << delta,
         0x1p+0;
    return B;
}

/// The FACTOR of the state weighting, not the weighting. Q = L^T L = I here,
/// and the visibility below is grouped so the cancellation happens against this
/// operand rather than against the assembled product.
auto family_state_weight_factor() -> weight_matrix
{
    weight_matrix L;
    L << 0x1p+0, 0x0p+0,
         0x0p+0, 0x1p+0;
    return L;
}

auto family_input_weight() -> input_weight
{
    input_weight R;
    R << 0x1p+0;
    return R;
}

auto enumerator_name(ctrlpp::dare_error error) -> const char *
{
    switch(error)
    {
        case ctrlpp::dare_error::non_stabilizable: return "non_stabilizable";
        case ctrlpp::dare_error::non_finite_input: return "non_finite_input";
        case ctrlpp::dare_error::singular_a:       return "singular_a";
        case ctrlpp::dare_error::singular_r:       return "singular_r";
        case ctrlpp::dare_error::singular_u11:     return "singular_u11";
        case ctrlpp::dare_error::non_psd_solution: return "non_psd_solution";
        case ctrlpp::dare_error::schur_failed:     return "schur_failed";
        case ctrlpp::dare_error::arithmetic_limit: return "arithmetic_limit";
    }
    return "unnamed";
}

/// Candidate quantity one: how visible the mode nearest the unit circle is
/// through the state weighting.
///
/// The discrete analogue of the continuous solver's nearest-imaginary-axis
/// visibility, with the unit circle in place of the axis. The Rayleigh quotient
/// v* Q v is grouped as ||L v||^2 against the weighting's factor L, because
/// evaluating it as a product against the assembled Q cancels the information
/// away where the quotient is many decades under the entries of Q and reports
/// exactly zero. Q real symmetric gives v* Q v = ||L Re v||^2 + ||L Im v||^2, so
/// no complex matrix product is needed. Dividing by the weighting's own
/// magnitude makes the result independent of how the caller scaled it.
auto unit_circle_mode_visibility(const state_matrix &A, const weight_matrix &weight_factor) -> double
{
    const Eigen::EigenSolver<state_matrix> eigensystem(A, true);

    double nearest       = std::numeric_limits<double>::infinity();
    int    nearest_index = 0;
    for(int index = 0; index < nx; ++index)
    {
        const double distance = std::abs(std::abs(eigensystem.eigenvalues()(index)) - 1.0);
        if(distance < nearest)
        {
            nearest       = distance;
            nearest_index = index;
        }
    }

    // `eigenvectors()` returns the matrix BY VALUE, so it must be named before a
    // column of it is taken. Writing `eigensystem.eigenvectors().col(index)`
    // deduces a block expression holding a reference into a temporary that dies
    // at the end of that statement, and every read below it is then undefined.
    // It does not fail loudly, and it is optimization-dependent, which is
    // precisely the axis this program sweeps.
    const Eigen::EigenSolver<state_matrix>::EigenvectorsType modes = eigensystem.eigenvectors();

    const auto          mode      = modes.col(nearest_index);
    const column_vector real_part = mode.real();
    const column_vector imag_part = mode.imag();

    const double quotient = (weight_factor * real_part).squaredNorm() + (weight_factor * imag_part).squaredNorm();

    // The denominator carries no cancellation -- it is a norm of nonnegative
    // contributions -- so it is taken from the assembled weighting.
    const double weight_magnitude = (weight_factor.transpose() * weight_factor).norm();
    return quotient / weight_magnitude;
}

/// The coupling of the mode outside the unit circle into the input.
///
/// The left eigenvector of A belonging to the eigenvalue of largest modulus,
/// applied to B. For this family that eigenvector is the first unit axis and the
/// coupling is |delta|, but that is a property of A and is measured here rather
/// than assumed.
auto outside_circle_mode_controllability(const state_matrix &A, const input_matrix &B) -> double
{
    const Eigen::EigenSolver<state_matrix> eigensystem(A.transpose(), true);

    double farthest      = -1.0;
    int    farthest_index = 0;
    for(int index = 0; index < nx; ++index)
    {
        const double modulus = std::abs(eigensystem.eigenvalues()(index));
        if(modulus > farthest)
        {
            farthest       = modulus;
            farthest_index = index;
        }
    }

    const Eigen::EigenSolver<state_matrix>::EigenvectorsType modes = eigensystem.eigenvectors();

    const auto          mode      = modes.col(farthest_index);
    const column_vector real_part = mode.real();
    const column_vector imag_part = mode.imag();

    const double real_coupling = (real_part.transpose() * B).norm();
    const double imag_coupling = (imag_part.transpose() * B).norm();
    return std::hypot(real_coupling, imag_coupling);
}

/// The symplectic matrix of the family, per Laub 1979 Eq. 7:
///
///     Z = [[A + G A^-T Q,  -G A^-T],
///          [   -A^-T Q,      A^-T ]]      with G = B R^-1 B^T
auto build_symplectic(const state_matrix  &A,
                      const input_matrix  &B,
                      const weight_matrix &Q,
                      const input_weight  &R,
                      symplectic_matrix   &Z) -> bool
{
    const state_matrix identity = state_matrix::Identity();

    const Eigen::ColPivHouseholderQR<transposed_state_matrix> qr_A_transpose = A.transpose().colPivHouseholderQr();
    if(!qr_A_transpose.isInvertible())
        return false;
    const state_matrix A_inverse_transpose = qr_A_transpose.solve(identity);

    const Eigen::ColPivHouseholderQR<input_weight> qr_R = R.colPivHouseholderQr();
    if(!qr_R.isInvertible())
        return false;
    const Eigen::Matrix<double, nu, nx> R_inverse_B_transpose = qr_R.solve(B.transpose());
    const state_matrix                  G                     = B * R_inverse_B_transpose;

    Z.block<nx, nx>(0, 0)   = A + G * A_inverse_transpose * Q;
    Z.block<nx, nx>(0, nx)  = -G * A_inverse_transpose;
    Z.block<nx, nx>(nx, 0)  = -A_inverse_transpose * Q;
    Z.block<nx, nx>(nx, nx) = A_inverse_transpose;

    return Z.allFinite();
}

/// Candidate quantity two: the rank-revealing pivot ratio of the leading block
/// of an orthonormal basis of the stable invariant subspace.
///
/// This is the quantity the solver's extraction compares against its rank
/// threshold, and it is reconstructed here by a route the solver does not take.
/// The basis is obtained by inverse orthogonal iteration rather than from
/// eigenvectors: the two stable eigenvalues of Z collide as the swept parameter
/// shrinks, and an eigenvector basis of a colliding pair loses half the
/// significand, while orthogonal iteration never forms one. Each step applies
/// Z^-1, which amplifies the stable directions relative to the rest by the ratio
/// of the two spectral radii, and reorthonormalizes.
///
/// The step count is not a tolerance and nothing is compared against it. It is
/// large enough that the iteration is stationary long before it ends, which the
/// report checks by printing the ratio at two different step counts.
auto stable_subspace_pivot_ratio(const symplectic_matrix &Z, int steps, double &ratio) -> bool
{
    const Eigen::FullPivLU<symplectic_matrix> factorization(Z);
    if(!factorization.isInvertible())
        return false;

    subspace_basis basis = subspace_basis::Zero();
    for(int index = 0; index < nx; ++index)
    {
        basis(index, index)      = 0x1p+0;
        basis(nx + index, index) = 0x1p+0;
    }

    for(int step = 0; step < steps; ++step)
    {
        const subspace_basis applied = factorization.solve(basis);
        if(!applied.allFinite())
            return false;

        const Eigen::HouseholderQR<subspace_basis> qr(applied);
        const symplectic_matrix                    q_factor = qr.householderQ();
        basis                                               = q_factor.leftCols(nx);
    }

    const state_matrix leading_block = basis.topRows(nx);

    // The solver factors the TRANSPOSE of this block, because it solves against
    // U11^T rather than inverting U11. The same operand is factored here.
    const Eigen::ColPivHouseholderQR<transposed_state_matrix> qr_leading = leading_block.transpose().colPivHouseholderQr();

    double smallest = std::numeric_limits<double>::infinity();
    double largest  = 0.0;
    for(int index = 0; index < nx; ++index)
    {
        const double pivot = std::abs(qr_leading.matrixQR()(index, index));
        if(pivot < smallest)
            smallest = pivot;
        if(pivot > largest)
            largest = pivot;
    }

    ratio = largest > 0.0 ? smallest / largest : 0.0;
    return true;
}

struct probe_point
{
    double             delta                 = 0.0;
    bool               accepted              = false;
    ctrlpp::dare_error error                 = ctrlpp::dare_error::non_stabilizable;
    double             visibility            = 0.0;
    double             pivot_ratio           = 0.0;
    bool               pivot_ratio_formed    = false;
    double             pivot_ratio_half      = 0.0;
    double             controllability       = 0.0;
    double             solution_magnitude    = 0.0;
    double             subspace_separation   = 0.0;
};

auto probe(double delta) -> probe_point
{
    probe_point point;
    point.delta = delta;

    const state_matrix  A = family_state();
    const input_matrix  B = family_input(delta);
    const weight_matrix L = family_state_weight_factor();
    const weight_matrix Q = (L.transpose() * L).eval();
    const input_weight  R = family_input_weight();

    const auto result = ctrlpp::dare<double, nx, nu>(A, B, Q, R);
    point.accepted    = result.has_value();
    if(result.has_value())
    {
        point.solution_magnitude  = result->P.cwiseAbs().maxCoeff();
        point.subspace_separation = result->subspace_separation;
    }
    else
    {
        point.error = result.error();
    }

    point.visibility      = unit_circle_mode_visibility(A, L);
    point.controllability = outside_circle_mode_controllability(A, B);

    symplectic_matrix Z = symplectic_matrix::Zero();
    if(build_symplectic(A, B, Q, R, Z))
    {
        double converged = 0.0;
        double halfway   = 0.0;
        if(stable_subspace_pivot_ratio(Z, 240, converged) && stable_subspace_pivot_ratio(Z, 120, halfway))
        {
            point.pivot_ratio        = converged;
            point.pivot_ratio_half   = halfway;
            point.pivot_ratio_formed = true;
        }
    }

    return point;
}

auto solve_family(double delta) -> ctrlpp::expected<ctrlpp::dare_result<double, nx, nu>, ctrlpp::dare_error>
{
    const state_matrix  A = family_state();
    const input_matrix  B = family_input(delta);
    const weight_matrix L = family_state_weight_factor();
    const weight_matrix Q = (L.transpose() * L).eval();
    const input_weight  R = family_input_weight();

    return ctrlpp::dare<double, nx, nu>(A, B, Q, R);
}

/// Predicate for the outer transition: the solve returned a solution.
auto accepts(double delta) -> bool
{
    const auto result = solve_family(delta);
    return result.has_value();
}

/// Predicate for the inner transition: the refusal is not the extraction's
/// rank verdict. True on the wide-delta side of that boundary, so the same
/// bisection drives both.
auto refusal_precedes_rank_verdict(double delta) -> bool
{
    const auto result = solve_family(delta);
    return result.has_value() || result.error() != ctrlpp::dare_error::singular_u11;
}

auto decade_exponent(double value) -> double
{
    return value > 0.0 ? std::log10(value) : -std::numeric_limits<double>::infinity();
}

void print_point(const char *label, const probe_point &point)
{
    std::printf("%-10s delta=%-22a log10_delta=%+9.4f status=%-7s enumerator=%-17s "
                "log10_visibility=%+9.4f log10_pivot_ratio=%+9.4f log10_pivot_ratio_half=%+9.4f "
                "log10_controllability=%+9.4f log10_solution_magnitude=%+9.4f subspace_separation=%.6e\n",
                label,
                point.delta,
                decade_exponent(point.delta),
                point.accepted ? "success" : "refused",
                point.accepted ? "-" : enumerator_name(point.error),
                decade_exponent(point.visibility),
                point.pivot_ratio_formed ? decade_exponent(point.pivot_ratio) : std::numeric_limits<double>::quiet_NaN(),
                point.pivot_ratio_formed ? decade_exponent(point.pivot_ratio_half) : std::numeric_limits<double>::quiet_NaN(),
                decade_exponent(point.controllability),
                point.accepted ? decade_exponent(point.solution_magnitude) : std::numeric_limits<double>::quiet_NaN(),
                point.accepted ? point.subspace_separation : std::numeric_limits<double>::quiet_NaN());
}

/// The decade ladder, pinned as hexadecimal literals so every compiler sweeps
/// the identical set of inputs. Decades rather than powers of two because the
/// enumerator's own documentation records its measured transition in decades.
constexpr double decade_ladder[] = {
    0x1.0000000000000p+0,  0x1.999999999999ap-4,  0x1.47ae147ae147bp-7,  0x1.0624dd2f1a9fcp-10,
    0x1.a36e2eb1c432dp-14, 0x1.4f8b588e368f1p-17, 0x1.0c6f7a0b5ed8dp-20, 0x1.ad7f29abcaf48p-24,
    0x1.5798ee2308c3ap-27, 0x1.12e0be826d695p-30, 0x1.b7cdfd9d7bdbbp-34, 0x1.5fd7fe1796495p-37,
    0x1.19799812dea11p-40, 0x1.c25c268497682p-44, 0x1.6849b86a12b9bp-47, 0x1.203af9ee75616p-50,
    0x1.cd2b297d889bcp-54,
};

/// The inputs the toolchain table is assembled over: bit-identical values, every
/// one of them a hexadecimal literal. Every configuration in the driver's matrix
/// evaluates these same literals, which is what makes the table a statement about
/// instruction selection rather than about the input.
///
/// The first four are the two crossings the bisection below resolves, each given
/// as the adjacent PAIR of doubles the verdict changes between. Anchoring the
/// table one representable number either side of a crossing is the strongest
/// placement available: it is where a difference in instruction selection has the
/// least distance to cover before it changes the answer. `deep` sits fifteen
/// binades below both, where the pivot ratio the rank test compares has reached
/// the arithmetic's own noise floor and has stopped falling with the input.
///
/// These four literals are the values one configuration's bisection returned.
/// They are pinned rather than recomputed precisely so that a configuration whose
/// crossing sits elsewhere reports a DIFFERENT verdict here, instead of quietly
/// moving the input to keep the verdict.
constexpr double accept_side_delta = 0x1.8801064100b01p-12;
constexpr double refuse_side_delta = 0x1.8801064100bp-12;
constexpr double accuracy_side_delta = 0x1.ed20087a00001p-25;
constexpr double rank_side_delta     = 0x1.ed20087ap-25;
constexpr double deep_anchor_delta   = 0x1p-40;

/// The state weighting shared by every input in this program, assembled from the
/// factor above rather than written out, so the factor stays the single source
/// the visibility is computed against.
auto family_state_weight() -> weight_matrix
{
    const weight_matrix L = family_state_weight_factor();
    return (L.transpose() * L).eval();
}

/// The second input the discrete error tests carry a disjunction over: a state
/// matrix and an input matrix that are both exactly zero. It is not a member of
/// the swept family, so it is solved on its own.
auto solve_zero_dynamics() -> ctrlpp::expected<ctrlpp::dare_result<double, nx, nu>, ctrlpp::dare_error>
{
    const state_matrix  A = state_matrix::Zero();
    const input_matrix  B = input_matrix::Zero();
    const weight_matrix Q = family_state_weight();
    const input_weight  R = family_input_weight();

    return ctrlpp::dare<double, nx, nu>(A, B, Q, R);
}

auto ulp_bits(double value) -> std::uint64_t
{
    return std::bit_cast<std::uint64_t>(value);
}

auto from_ulp_bits(std::uint64_t bits) -> double
{
    return std::bit_cast<double>(bits);
}

auto step_ulps(double value, int offset) -> double
{
    const std::uint64_t bits = ulp_bits(value);
    return from_ulp_bits(offset >= 0 ? bits + static_cast<std::uint64_t>(offset)
                                     : bits - static_cast<std::uint64_t>(-offset));
}

void print_neighborhood(const char *label, double center, int radius)
{
    for(int offset = -radius; offset <= radius; ++offset)
    {
        const probe_point point = probe(step_ulps(center, offset));
        std::printf("%-10s offset=%+3d ", label, offset);
        print_point("", point);
    }
}

auto verdict_word(const ctrlpp::expected<ctrlpp::dare_result<double, nx, nu>, ctrlpp::dare_error> &result) -> const char *
{
    return result.has_value() ? "success" : "refused";
}

auto verdict_enumerator(const ctrlpp::expected<ctrlpp::dare_result<double, nx, nu>, ctrlpp::dare_error> &result) -> const char *
{
    return result.has_value() ? "-" : enumerator_name(result.error());
}

/// Find the adjacent pair of ladder rungs a predicate changes across.
auto bracket_on_binary_ladder(bool (*holds)(double), int first_exponent, int last_exponent, double &upper, double &lower) -> bool
{
    upper = 0.0;
    lower = 0.0;
    for(int exponent = first_exponent; exponent <= last_exponent; ++exponent)
    {
        const double delta = std::ldexp(0x1p+0, -exponent);
        if(holds(delta))
        {
            if(lower == 0.0)
                upper = delta;
        }
        else if(upper > 0.0 && lower == 0.0)
        {
            lower = delta;
        }
    }
    return upper > 0.0 && lower > 0.0;
}

/// Bisect a predicate that holds at `upper` and fails at `lower` until the two
/// are adjacent doubles. The IEEE bit pattern of a positive double is monotone
/// in its value, so the crossing resolves to a single representable number.
///
/// Nothing about monotonicity is assumed beyond the two endpoints. The census
/// printed alongside is what reports whether the verdict changes once across the
/// neighborhood or interleaves inside it.
void bisect_crossing(double upper, double lower, bool (*holds)(double), double &last_true, double &first_false)
{
    std::uint64_t high = ulp_bits(upper);
    std::uint64_t low  = ulp_bits(lower);
    while(high - low > 1)
    {
        const std::uint64_t middle = low + (high - low) / 2;
        if(holds(from_ulp_bits(middle)))
            high = middle;
        else
            low = middle;
    }
    last_true   = from_ulp_bits(high);
    first_false = from_ulp_bits(low);
}

/// Locate a predicate's crossing from the binary ladder alone, with no pinned
/// bracket, so each configuration reports the crossing IT produces.
auto locate_crossing(bool (*holds)(double), double &last_true, double &first_false) -> bool
{
    double upper = 0.0;
    double lower = 0.0;
    if(!bracket_on_binary_ladder(holds, 4, 40, upper, lower))
        return false;
    bisect_crossing(upper, lower, holds, last_true, first_false);
    return true;
}

/// How far a crossing sits from a reference crossing, in decades. Reported this
/// way rather than in representable doubles because the two are hundreds of
/// millions of millions of doubles apart, a count that carries no meaning; the
/// shift is a fraction of a decade and that is the quantity worth reading.
auto shift_decades(double value, double reference) -> double
{
    if(!(value > 0.0) || !(reference > 0.0))
        return std::numeric_limits<double>::quiet_NaN();
    return std::log10(value) - std::log10(reference);
}

/// One line carrying every input the driver builds its table from, so a whole
/// configuration is one invocation and one row.
///
/// The seven verdict fields are evaluated at the pinned literals above. The two
/// crossing fields are each configuration's OWN bisected crossing, reported both
/// in hexadecimal and as its signed distance in decades from the pinned literal.
/// A configuration whose crossing sits exactly where the pinned one does reports
/// a shift of zero.
void report_anchor_line()
{
    const auto accept_side   = solve_family(accept_side_delta);
    const auto refuse_side   = solve_family(refuse_side_delta);
    const auto accuracy_side = solve_family(accuracy_side_delta);
    const auto rank_side     = solve_family(rank_side_delta);
    const auto deep_result   = solve_family(deep_anchor_delta);
    const auto family_zero   = solve_family(0x0p+0);
    const auto zero_dynamics = solve_zero_dynamics();

    double accept_last  = 0.0;
    double accept_first = 0.0;
    const bool accept_located = locate_crossing(&accepts, accept_last, accept_first);

    double rank_last  = 0.0;
    double rank_first = 0.0;
    const bool rank_located = locate_crossing(&refusal_precedes_rank_verdict, rank_last, rank_first);

    std::printf("accept_side=%s:%s refuse_side=%s:%s accuracy_side=%s:%s rank_side=%s:%s "
                "deep=%s:%s family_zero=%s:%s zero_dynamics=%s:%s "
                "accept_crossing=%a accept_shift=%+.6f "
                "rank_crossing=%a rank_shift=%+.6f\n",
                verdict_word(accept_side),
                verdict_enumerator(accept_side),
                verdict_word(refuse_side),
                verdict_enumerator(refuse_side),
                verdict_word(accuracy_side),
                verdict_enumerator(accuracy_side),
                verdict_word(rank_side),
                verdict_enumerator(rank_side),
                verdict_word(deep_result),
                verdict_enumerator(deep_result),
                verdict_word(family_zero),
                verdict_enumerator(family_zero),
                verdict_word(zero_dynamics),
                verdict_enumerator(zero_dynamics),
                accept_first,
                accept_located ? shift_decades(accept_first, refuse_side_delta) : std::numeric_limits<double>::quiet_NaN(),
                rank_first,
                rank_located ? shift_decades(rank_first, rank_side_delta) : std::numeric_limits<double>::quiet_NaN());
}

void report_transition(const char *name, double upper, double lower, bool (*holds)(double), int radius)
{
    double last_true   = 0.0;
    double first_false = 0.0;
    bisect_crossing(upper, lower, holds, last_true, first_false);

    std::printf("\n== %s: bisected to a single representable input ==\n", name);
    print_point("above", probe(last_true));
    print_point("below", probe(first_false));

    int held_below   = 0;
    int failed_below = 0;
    int held_above   = 0;
    int failed_above = 0;
    for(int offset = -radius; offset <= radius; ++offset)
    {
        if(offset == 0)
            continue;
        const bool value = holds(step_ulps(last_true, offset));
        if(offset < 0)
            (value ? held_below : failed_below) += 1;
        else
            (value ? held_above : failed_above) += 1;
    }
    std::printf("%s neighborhood radius=%d: smaller side held=%d failed=%d, larger side held=%d failed=%d\n",
                name,
                radius,
                held_below,
                failed_below,
                held_above,
                failed_above);
    std::printf("%s crossing: holds at %a (log10 %+.6f), fails at the adjacent double %a (log10 %+.6f)\n",
                name,
                last_true,
                decade_exponent(last_true),
                first_false,
                decade_exponent(first_false));
}

void report_full()
{
    std::printf("== decade ladder over the recorded family ==\n");
    for(const double delta : decade_ladder)
    {
        const probe_point point = probe(delta);
        print_point("decade", point);
    }

    std::printf("\n== binary ladder spanning both transitions ==\n");
    for(int exponent = 4; exponent <= 40; ++exponent)
    {
        const probe_point point = probe(std::ldexp(0x1p+0, -exponent));
        print_point("binary", point);
    }

    std::printf("\n== plus and minus three units in the last place at four decades ==\n");
    print_neighborhood("ulp-1e-3", 0x1.0624dd2f1a9fcp-10, 3);
    print_neighborhood("ulp-1e-4", 0x1.a36e2eb1c432dp-14, 3);
    print_neighborhood("ulp-1e-7", 0x1.ad7f29abcaf48p-24, 3);
    print_neighborhood("ulp-1e-8", 0x1.5798ee2308c3ap-27, 3);

    double upper = 0.0;
    double lower = 0.0;
    if(bracket_on_binary_ladder(&accepts, 4, 40, upper, lower))
        report_transition("accept-refuse", upper, lower, &accepts, 16);
    else
        std::printf("\naccept-refuse crossing: NOT bracketed inside the swept range\n");

    if(bracket_on_binary_ladder(&refusal_precedes_rank_verdict, 4, 40, upper, lower))
        report_transition("rank-verdict", upper, lower, &refusal_precedes_rank_verdict, 16);
    else
        std::printf("\nrank-verdict crossing: NOT bracketed inside the swept range\n");

    std::printf("\n== the two inputs the discrete error tests already carry a disjunction over ==\n");
    print_point("family-0", probe(0x0p+0));
    const auto zero_dynamics = solve_zero_dynamics();
    std::printf("%-10s status=%-7s enumerator=%s\n", "zero-dyn", verdict_word(zero_dynamics), verdict_enumerator(zero_dynamics));

    std::printf("\n== toolchain anchors ==\n");
    report_anchor_line();
}

}

auto main(int argc, char **argv) -> int
{
    if(argc > 1 && std::strcmp(argv[1], "--anchor") == 0)
    {
        report_anchor_line();
        return 0;
    }

    report_full();
    return 0;
}
