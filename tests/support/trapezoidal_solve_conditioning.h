#ifndef HPP_GUARD_CTRLPP_TESTS_SUPPORT_TRAPEZOIDAL_SOLVE_CONDITIONING_H
#define HPP_GUARD_CTRLPP_TESTS_SUPPORT_TRAPEZOIDAL_SOLVE_CONDITIONING_H

/// @brief One conditioning model for the trapezoidal retiming solve, shared by
/// the fuzz target and the anchor test.
///
/// It lived in two copies before, one per consumer, and the two had already
/// drifted apart on the answer they gave when the quantity they divide by was
/// exactly zero. A model that exists twice is a model that can be repaired once,
/// and the claim resting on it then rests on whichever copy was not touched.
/// There is one copy now and both consumers include it.

#include <cmath>
#include <limits>
#include <algorithm>

namespace ctrlpp::test
{

/// @brief The answer when an amplification has no finite value.
///
/// Never an allowance. Both consumers turn this into a failure: an accepted
/// solve whose amplification is unbounded is a solve whose result carries no
/// significant digits, which means the library returned success on a request it
/// could not resolve.
inline constexpr double unbounded_solve_conditioning = std::numeric_limits<double>::infinity();

/// @brief Which of the five expressions the retiming solve actually evaluated.
///
/// Three shapes, and two of them carry two forms: the plateau falls back to the
/// cruise velocity itself when its boundary velocity is not positive and there
/// is no boundary to sit near, and the valley carries whichever of the
/// cruise-velocity decrement and the cruise velocity is the smaller.
enum class trapezoidal_solve_form
{
    plateau_rise,
    plateau_cruise_velocity,
    ramp_through,
    valley_decrement,
    valley_cruise_velocity,
};

/// @brief Chained rounding operations behind the solved cruise velocity.
///
/// RE-DERIVED for the shifted parametrization; the count this replaced was
/// twelve and described a discriminant in the cruise velocity that neither
/// solved branch forms any more. The longest of the five paths above governs,
/// because a budget has to bound from above, and that is the shifted one:
/// fifteen operations form the boundary duration the shift is measured from
/// (the sibling count the library carries for its own resolution floor), six
/// form the ramp-through residual (five for the swept ramp distance and a sixth
/// for the subtraction), one forms the duration increment, three form the
/// linear coefficient (the residual over the boundary velocity, the sum, and
/// the multiplication by the acceleration), two form the constant term, three
/// form the discriminant, one takes its square root, three form the
/// addition-only root selection (the sum in the denominator, the doubled
/// constant, and the division), and one forms the cruise velocity from the
/// boundary velocity and the root. Fifteen plus six plus one plus three plus
/// two plus three plus one plus three plus one is thirty-five.
///
/// The other four paths are shorter and are bounded by it: the ramp-through
/// shape spends ten, and both cruise-velocity forms spend sixteen on the value
/// they return. Every operation is counted whether or not it rounds, so the
/// count bounds the accumulated error from above rather than describing it
/// tightly.
inline constexpr int trapezoidal_cruise_velocity_rounding_ops = 35;

/// @brief Chained rounding operations behind the duration recomputed from that
/// velocity: a subtraction and a division for each of the two ramp durations
/// (four), an addition, a multiplication and a division for each of the two
/// swept ramp distances (six more, ten), two subtractions and a division for
/// the cruise term (three more, thirteen), and two additions for the final sum.
inline constexpr int trapezoidal_realized_duration_rounding_ops = 15;

/// @brief The two chains above, which together stand behind a realized duration.
inline constexpr int trapezoidal_duration_rounding_ops =
    trapezoidal_cruise_velocity_rounding_ops + trapezoidal_realized_duration_rounding_ops;

/// @brief A total duration paired with the largest operand that entered it.
template <typename Scalar>
struct trapezoidal_boundary_duration
{
    Scalar T{};
    Scalar scale{};
};

/// @brief The duration at a cruise velocity and the scale its own rounding is
/// measured against, recomputed exactly as the library forms them.
///
/// The scale is not the duration alone: the cruise term's numerator is the
/// commanded displacement less the two swept ramp distances, and that difference
/// cancels to nothing whenever the ramps sweep almost all of the displacement.
template <typename Scalar>
auto trapezoidal_duration_and_scale_at(Scalar v, Scalar v0, Scalar v1, Scalar a, Scalar h)
    -> trapezoidal_boundary_duration<Scalar>
{
    Scalar const T_a = std::abs(v - v0) / a;
    Scalar const T_d = std::abs(v - v1) / a;
    Scalar const d_a = (v0 + v) * T_a / Scalar{2};
    Scalar const d_d = (v1 + v) * T_d / Scalar{2};
    Scalar const T = T_a + (h - d_a - d_d) / v + T_d;
    Scalar const cruise_scale = std::max({h, std::abs(d_a), std::abs(d_d)}) / std::abs(v);
    return {T, std::max({T_a, T_d, cruise_scale, std::abs(T)})};
}

/// @brief Relative-error amplification of a difference: the scale of the
/// operands that entered it over the magnitude of what came out.
///
/// A result of exactly zero has no relative error to speak of and the honest
/// answer is that the amplification is unbounded. Neither of the two answers the
/// retired copies gave here -- zero in one, one in the other -- is that.
inline auto trapezoidal_amplification(double operands, double result) -> double
{
    double const magnitude = std::abs(result);
    return (magnitude > 0.0) ? std::abs(operands) / magnitude : unbounded_solve_conditioning;
}

/// @brief Absolute error of the square root of a quadratic's discriminant, in
/// units of the scalar type's epsilon.
///
/// The two coefficients arrive as ABSOLUTE error budgets -- the scale of the
/// operands that entered each -- rather than as the relative amplifications this
/// took before. The distinction is not presentational. Every term below needs
/// the absolute budget, so a relative form has to be multiplied back by the very
/// magnitude it was divided by, and that product is zero times infinity exactly
/// where a coefficient cancels completely. A total cancellation is the case this
/// model exists to describe, not one it may answer NaN on.
///
/// The magnitudes are taken as magnitudes for the same reason. An error bound
/// built from a signed coefficient is reduced by that coefficient's sign, and
/// the linear one is negative over half the domain of the root selection below.
///
/// The discriminant is `b^2 - 4c`, so its own absolute error is
/// `2 |b| E_b + 4 E_c` with the two budgets given. Two bounds on the root's
/// error follow from it and both hold, so the smaller governs:
///
///  * the LINEARIZED one, half the discriminant's error over the root. It is the
///    right answer while that error is small against the discriminant.
///  * the MERGED one, the square root of the discriminant's own uncertainty.
///    Where a quadratic's two roots merge the discriminant is zero to the
///    precision available, the linearization diverges, and this is what remains:
///    the root cannot be larger than the square root of the largest the
///    discriminant could have been.
///
/// The merged bound matters rather than being a formality. A discriminant that
/// rounds to exactly zero is a request sitting on its shape's reachable extreme,
/// where the root selection returns the doubled constant over the linear
/// coefficient with nothing to cancel; the linearized answer would call that
/// unbounded and, since an unbounded answer is a failure, would report a defect
/// where the solve carries one unit in the last place.
///
/// Both are formed by dividing before multiplying, and the merged one by
/// splitting the square root across the two terms it sums -- the square root of
/// a sum is at most the sum of the square roots, and its linear term is split
/// once more across its own two factors. The discriminant's error is never
/// FORMED, only ever divided into: it is the product of a coefficient with a
/// budget carried at that coefficient's own scale, and it leaves the
/// representable range on requests the solve resolves without difficulty. A
/// long move at a small boundary velocity reaches it while every quantity the
/// model is built from is still ordinary.
///
/// A linearized term is left NaN rather than guarded where the root it divides
/// by is zero and the budget it divides is too. `fmin` returns the other operand
/// when one is NaN, so the merged bound governs there -- which is precisely the
/// case the merged bound exists for, and a guard would only spell it a second
/// time.
template <typename Scalar>
auto trapezoidal_discriminant_root_error(double linear_error, Scalar linear, double constant_error,
                                         Scalar discriminant) -> double
{
    constexpr double eps = static_cast<double>(std::numeric_limits<Scalar>::epsilon());
    double const b = std::abs(static_cast<double>(linear));
    double const root = std::sqrt(std::max(static_cast<double>(discriminant), 0.0));

    double const linearized =
        0.5 * (2.0 * (b / root) * linear_error + 4.0 * (constant_error / root));
    double const merged = std::sqrt(2.0 * linear_error / eps) * std::sqrt(b)
                          + std::sqrt(4.0 * (constant_error / eps));
    return std::fmin(linearized, merged);
}

/// @brief Every intermediate a shifted solve forms, recomputed in the scalar
/// type the library solved in.
///
/// The recomputation is in the profile's own scalar type on purpose. A model
/// evaluated in double precision against a single-precision solve sees a
/// different cancellation than the one that happened, and the difference is
/// worst exactly where the shift matters: the increment measured from a boundary
/// duration can be several units in the last place in the type that solved it
/// and round to zero in a wider one. Mirroring the type is what makes an
/// unbounded answer here mean something about the library rather than about the
/// model.
template <typename Scalar>
struct trapezoidal_shifted_solve
{
    Scalar dT{};    ///< duration increment (valley) or decrement (plateau)
    Scalar scale{}; ///< the boundary duration's own operand scale
    Scalar B{};
    Scalar C{};
    Scalar disc{};
    Scalar root{}; ///< the decrement below (valley) or the rise above (plateau)
    bool solved{};
};

template <typename Scalar>
auto trapezoidal_trace_shifted_solve(bool valley, Scalar a, Scalar h, Scalar v0, Scalar v1,
                                     Scalar T_target) -> trapezoidal_shifted_solve<Scalar>
{
    Scalar const v_lo = std::min(v0, v1);
    Scalar const v_hi = std::max(v0, v1);
    Scalar const v_ref = valley ? v_lo : v_hi;
    Scalar const ramp_residual = h - (v_hi - v_lo) * (v_hi + v_lo) / (Scalar{2} * a);
    auto const at = trapezoidal_duration_and_scale_at(v_ref, v0, v1, a, h);

    trapezoidal_shifted_solve<Scalar> trace;
    trace.scale = at.scale;
    trace.dT = valley ? (T_target - at.T) : (at.T - T_target);
    trace.B = valley ? a * (trace.dT + ramp_residual / v_lo)
                     : a * (ramp_residual / v_hi - trace.dT);
    trace.C = a * v_ref * trace.dT;
    trace.disc = trace.B * trace.B - Scalar{4} * trace.C;
    trace.solved = (trace.B > Scalar{0}) && (trace.disc >= Scalar{0});
    trace.root = trace.solved ? Scalar{2} * trace.C / (trace.B + std::sqrt(trace.disc)) : Scalar{0};
    return trace;
}

/// @brief Which expression the solve took, decided the way the library decides.
///
/// The shape test is the one the library uses to ACCEPT a solved cruise
/// velocity, which admits both boundary velocities, rather than the
/// strict-interior test the retired copies used. The two disagree at the smaller
/// boundary velocity, and there the retired model applied the valley formula to
/// a ramp-through solve.
///
/// The valley's two forms are separated by the library's own comparison, which
/// asks only which of the decrement and the cruise velocity it leaves behind is
/// larger and introduces no constant.
template <typename Scalar>
auto trapezoidal_solve_form_taken(Scalar v_cruise, Scalar a, Scalar h, Scalar v0, Scalar v1,
                                  Scalar T_target) -> trapezoidal_solve_form
{
    Scalar const v_lo = std::min(v0, v1);
    Scalar const v_hi = std::max(v0, v1);
    if(v_cruise >= v_hi)
        return (v_hi > Scalar{0}) ? trapezoidal_solve_form::plateau_rise
                                  : trapezoidal_solve_form::plateau_cruise_velocity;
    if(v_cruise >= v_lo)
        return trapezoidal_solve_form::ramp_through;

    auto const trace = trapezoidal_trace_shifted_solve(true, a, h, v0, v1, T_target);
    return (trace.solved && trace.root + trace.root <= v_lo)
               ? trapezoidal_solve_form::valley_decrement
               : trapezoidal_solve_form::valley_cruise_velocity;
}

/// @brief Sensitivity of the realized duration to the cruise velocity, as a
/// relative-to-relative ratio: |dT/dv| * v / T.
///
/// Within a shape T(v) = (+-)(v - v0 - v1)/a + D/v with D that shape's residual
/// displacement, so |dT/dv| * v is bounded by v/a + |D|/v. On the plateau that
/// bound is loose by orders of magnitude, because the derivative vanishes at the
/// triangular cruise velocity and the two terms very nearly cancel there; the
/// exact magnitude |v_tri - v|(v_tri + v)/(a v) is used instead, and it is safe
/// to form because the triangular velocity is a square root of a SUM of
/// nonnegative terms. The valley has no such stationary point available: its
/// counterpart is the vanishing-cruise velocity, a square root of a difference
/// of two nearly equal quantities that carries no significant digits in exactly
/// the configuration the valley solve has to get right, so the valley keeps the
/// bound.
template <typename Scalar>
auto trapezoidal_duration_sensitivity(trapezoidal_solve_form form, Scalar v_cruise, Scalar a,
                                      Scalar h, Scalar v0, Scalar v1, Scalar T_target) -> double
{
    Scalar const v_lo = std::min(v0, v1);
    Scalar const v_hi = std::max(v0, v1);
    Scalar const v_sum_sq = v0 * v0 + v1 * v1;
    double const v = std::abs(static_cast<double>(v_cruise));
    double slope_times_v = 0.0;

    switch(form)
    {
    case trapezoidal_solve_form::plateau_rise:
    case trapezoidal_solve_form::plateau_cruise_velocity:
    {
        auto const v_tri = static_cast<double>(std::sqrt(a * h + v_sum_sq / Scalar{2}));
        slope_times_v = std::abs(v_tri - v) * (v_tri + v) / (static_cast<double>(a) * v);
        break;
    }
    case trapezoidal_solve_form::ramp_through:
    {
        // Both ramps run one way, so their duration and their swept distance are
        // both independent of the cruise velocity and only the residual remains.
        Scalar const ramp_residual = h - (v_hi - v_lo) * (v_hi + v_lo) / (Scalar{2} * a);
        slope_times_v = std::abs(static_cast<double>(ramp_residual)) / v;
        break;
    }
    default:
    {
        double const residual_displacement =
            static_cast<double>(h) - static_cast<double>(v_sum_sq) / (2.0 * static_cast<double>(a));
        slope_times_v = v / static_cast<double>(a) + std::abs(residual_displacement) / v;
        break;
    }
    }

    return (T_target > Scalar{0}) ? slope_times_v / static_cast<double>(T_target)
                                  : unbounded_solve_conditioning;
}

/// @brief Relative-error amplification of the solved cruise velocity itself, in
/// units of the scalar type's epsilon.
///
/// Each form reaches its answer through at least one difference of two nearly
/// equal quantities, and the ratio of the operands entering that difference to
/// the difference itself is what the answer inherits. The differences are named
/// per form:
///
///  * both shifted forms: the ramp-through residual (the commanded displacement
///    less the distance the two ramps sweep on their own), the duration
///    increment or decrement measured from a boundary duration, the linear
///    coefficient built from those two, and the discriminant. The final step
///    adds or subtracts the root from the boundary velocity, and its
///    contribution is the root's own relative error scaled by how small the root
///    is against the velocity it sits on -- which is the whole point of solving
///    for the root rather than the velocity.
///  * ramp-through: the residual again, and the duration less the fixed duration
///    of the two ramps.
///  * both cruise-velocity forms: the discriminant, and -- only where the linear
///    coefficient is positive and the root selection therefore divides by it --
///    the constant term, which is the residual displacement carried at the scale
///    of the squared boundary velocities. Where that coefficient is not positive
///    the library adds two magnitudes instead of dividing, and the constant term
///    contributes through the discriminant alone.
template <typename Scalar>
auto trapezoidal_velocity_amplification(trapezoidal_solve_form form, Scalar v_cruise, Scalar a,
                                        Scalar h, Scalar v0, Scalar v1, Scalar T_target) -> double
{
    Scalar const v_lo = std::min(v0, v1);
    Scalar const v_hi = std::max(v0, v1);
    Scalar const v_sum_sq = v0 * v0 + v1 * v1;
    double const v = std::abs(static_cast<double>(v_cruise));

    // The swept ramp distance is formed factored, as a difference of the two
    // boundary velocities times their sum, and one of those two is exact
    // whenever the other cancels: two values within a factor of two of each
    // other subtract exactly, and outside that range there is nothing to cancel.
    // So the distance carries a handful of units in the last place at its OWN
    // scale, and the operands of the residual are the displacement and the
    // distance rather than the squared velocities the library's own rejection
    // floor is conservatively written against.
    Scalar const ramp_distance = (v_hi - v_lo) * (v_hi + v_lo) / (Scalar{2} * a);
    Scalar const ramp_residual = h - ramp_distance;
    double const residual_error = std::max(std::abs(static_cast<double>(h)),
                                           std::abs(static_cast<double>(ramp_distance)));

    switch(form)
    {
    case trapezoidal_solve_form::plateau_rise:
    case trapezoidal_solve_form::valley_decrement:
    {
        bool const valley = (form == trapezoidal_solve_form::valley_decrement);
        Scalar const v_ref = valley ? v_lo : v_hi;
        auto const trace = trapezoidal_trace_shifted_solve(valley, a, h, v0, v1, T_target);

        // The residual reaches the linear coefficient divided by the boundary
        // velocity, so its budget is carried through that same division. A
        // boundary velocity of exactly zero leaves it unbounded, which is the
        // honest answer and is unreachable on an accepted retiming: the cruise
        // velocity a shifted solve returns is that boundary offset by the root,
        // and a retiming whose cruise velocity is not strictly positive is
        // rejected before it is returned.
        double const increment_error = std::abs(static_cast<double>(trace.scale));
        double const boundary = std::abs(static_cast<double>(v_ref));
        double const coefficient_error =
            static_cast<double>(a)
            * (trapezoidal_amplification(residual_error, v_ref) + increment_error);
        double const constant_error = static_cast<double>(a) * boundary * increment_error;
        // The root selection ADDS the square root to the linear coefficient
        // rather than subtracting it, so its denominator carries the two
        // absolute errors side by side and neither is amplified by the other's
        // smallness.
        double const coefficient = std::abs(static_cast<double>(trace.B));
        double const denominator_error =
            coefficient_error
            + trapezoidal_discriminant_root_error(coefficient_error, trace.B, constant_error,
                                                  trace.disc);
        double const denominator =
            coefficient + std::sqrt(std::max(static_cast<double>(trace.disc), 0.0));
        // The root's ABSOLUTE error, formed without ever passing through its
        // relative one. The constant term's budget reaches it doubled over the
        // same denominator the root itself carries, because the root IS the
        // doubled constant term over that denominator; the denominator's budget
        // reaches it scaled by the root. Dividing by the cruise velocity at the
        // end turns it back into the relative amplification this function
        // reports, and that division is the only one taken.
        double const root_error = (2.0 * constant_error
                                   + std::abs(static_cast<double>(trace.root)) * denominator_error)
                                  / denominator;
        return root_error / v;
    }
    case trapezoidal_solve_form::ramp_through:
    {
        Scalar const ramp_span = std::abs(v1 - v0) / a;
        Scalar const denominator = T_target - ramp_span;
        return trapezoidal_amplification(residual_error, ramp_residual)
               + trapezoidal_amplification(std::max(T_target, ramp_span), denominator);
    }
    default:
    {
        bool const plateau = (form == trapezoidal_solve_form::plateau_cruise_velocity);
        Scalar const b = plateau ? ((v0 + v1) + a * T_target) : (a * T_target - (v0 + v1));
        Scalar const c = plateau ? (a * h + v_sum_sq / Scalar{2}) : (v_sum_sq / Scalar{2} - a * h);
        double const linear_error =
            std::max(std::abs(static_cast<double>(a * T_target)),
                     std::abs(static_cast<double>(v0 + v1)));
        double const constant_error = std::max(std::abs(static_cast<double>(a * h)),
                                               std::abs(static_cast<double>(v_sum_sq / Scalar{2})));
        // Both spellings of the root selection add magnitudes rather than
        // subtracting them, so the same side-by-side treatment applies.
        Scalar const disc = b * b - Scalar{4} * c;
        double const linear = std::abs(static_cast<double>(b));
        double const denominator_error =
            linear_error
            + trapezoidal_discriminant_root_error(linear_error, b, constant_error, disc);
        double const denominator = linear + std::sqrt(std::max(static_cast<double>(disc), 0.0));
        double const selection = trapezoidal_amplification(denominator_error, denominator);
        // The library spells this root selection two ways and branches on the
        // SIGN of the linear coefficient, so the model has to branch with it. A
        // positive coefficient takes the rationalized spelling, which divides the
        // doubled constant term by the sum and inherits that term's relative
        // error. A non-positive one takes the direct spelling, which adds the
        // root to the coefficient's magnitude and halves the sum: two magnitudes
        // added, no division by the constant term, and none of its relative error
        // inherited. Charging that error to both spellings called the solve
        // unbounded wherever the constant term cancels exactly -- a displacement
        // the two ramps sweep precisely, which the shifted parametrization is
        // built to handle -- on retimings the library resolves to the last bit.
        return (b > Scalar{0}) ? trapezoidal_amplification(constant_error, c) + selection
                               : selection;
    }
    }
}

/// @brief Amplification the retiming applies to its own roundings, as a
/// multiplier on a budget of units in the last place at the requested duration.
///
/// The realized duration is recomputed from the cruise velocity the solve
/// returns, so two independent things can amplify a rounding into a duration
/// error and the larger of them governs.
///
/// The first is the solve: the velocity's relative error, times the duration's
/// sensitivity to it. The second is the recomputation itself, which forms the
/// cruise phase as the commanded displacement less the two swept ramp distances
/// divided by the velocity -- a difference that cancels to nothing whenever the
/// ramps sweep almost all of the displacement, and a division that magnifies
/// whatever survives when the cruise velocity is small. That second term is
/// intrinsic to the profile family and to no parametrization of the solve: a
/// cruise velocity many orders of magnitude below the boundary velocities
/// divides a residual by an almost-zero number, and no rearrangement of the
/// solve changes what that does. It is the same scale the library carries
/// alongside its own boundary durations, so the two cannot drift apart.
///
/// An amplification with no finite value is reported as such and is never turned
/// into an allowance. Every quantity the model divides by is one the library
/// itself floored before it accepted the solve, recomputed in the same scalar
/// type, so an unbounded answer on an accepted retiming says the library
/// returned success where it had no digits to return it with.
///
/// Returning the amplification rather than folding a fixed factor in keeps the
/// bound tight where the problem is well conditioned, which is where a defect
/// would actually have to hide.
template <typename Scalar>
auto trapezoidal_solve_conditioning(Scalar v_cruise, Scalar a, Scalar h, Scalar v0, Scalar v1,
                                    Scalar T_target) -> double
{
    auto const form = trapezoidal_solve_form_taken(v_cruise, a, h, v0, v1, T_target);
    double const velocity =
        trapezoidal_velocity_amplification(form, v_cruise, a, h, v0, v1, T_target);
    double const sensitivity =
        trapezoidal_duration_sensitivity(form, v_cruise, a, h, v0, v1, T_target);

    auto const realized = trapezoidal_duration_and_scale_at(v_cruise, v0, v1, a, h);
    double const recomputation = trapezoidal_amplification(realized.scale, T_target);

    if(!std::isfinite(velocity) || !std::isfinite(sensitivity) || !std::isfinite(recomputation))
        return unbounded_solve_conditioning;
    return std::max({1.0, velocity * sensitivity, recomputation});
}

}

#endif
