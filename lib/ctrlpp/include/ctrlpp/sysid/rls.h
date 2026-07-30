#ifndef HPP_GUARD_CTRLPP_SYSID_RLS_H
#define HPP_GUARD_CTRLPP_SYSID_RLS_H

/// @brief Recursive Least Squares with bounded covariance and forgetting factor.
///
/// @cite ljung1999 -- Ljung, "System Identification: Theory for the User", 1999

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <cstddef>
#include <utility>
#include <algorithm>

namespace ctrlpp
{

/// @brief Structured failure modes of the `rls` factory.
///
/// The first enumerator is an exact domain condition read off the covariance
/// update; the second enforces the estimator's documented parameter contract
/// and says so, rather than dressing a convention up as arithmetic.
///
///  * non_positive_forgetting_factor : the forgetting factor is not finite and
///                                     strictly positive. It is the DIVISOR of
///                                     the covariance update
///                                     P <- (P - k phi^T P) / lambda, so this is
///                                     read off the arithmetic and is not a
///                                     preference: at zero the update divides by
///                                     zero and every entry of P becomes
///                                     non-finite on the first sample, and at a
///                                     negative value the division negates a
///                                     positive semidefinite matrix, so every
///                                     later gain points against the error
///                                     instead of along it. A non-finite factor
///                                     poisons P by the same division.
///  * forgetting_factor_above_unity  : the forgetting factor exceeds one. The
///                                     arithmetic does NOT forbid this and the
///                                     distinction is worth stating: dividing by
///                                     a factor above one DEFLATES the
///                                     covariance faster than the measurement
///                                     update alone would, so the gain collapses
///                                     toward zero and the estimator silently
///                                     stops adapting -- it does not diverge.
///                                     Exponential forgetting is defined on
///                                     (0, 1], which is the contract this type
///                                     documents, and a value above one is
///                                     almost always an inverted reading of the
///                                     parameter. It is rejected as a contract
///                                     violation, not as a domain violation.
///  * non_finite_initial_covariance  : the initial covariance has a non-finite
///                                     entry. It seeds the recursion and the
///                                     recursion has no mechanism that could
///                                     return a non-finite covariance to a
///                                     finite one. Finiteness is the whole test:
///                                     a singular or zero P0 is a legitimate,
///                                     deliberately posed starting point and is
///                                     accepted.
///  * non_positive_covariance_bound  : the covariance upper bound is not finite
///                                     and strictly positive. When the trace
///                                     exceeds NP times the bound the update
///                                     rescales P by bound*NP/trace, so a
///                                     negative bound negates the covariance and
///                                     a zero bound drives it to exactly zero on
///                                     the first sample, after which the gain is
///                                     zero forever and the parameters never
///                                     move again. A non-finite bound disables
///                                     the clamp the field exists to impose.
enum class rls_error
{
    non_positive_forgetting_factor,
    forgetting_factor_above_unity,
    non_finite_initial_covariance,
    non_positive_covariance_bound,
};

/// @brief Structured refusal modes of a single `rls::update` cycle.
///
/// The order below is the order the guard tests them, and it is a severity
/// order: the estimator's own carried state first, then the caller's two
/// operands, then what the arithmetic formed from them. The rule is that the
/// most upstream cause is named, because a caller told "your regressor carries
/// no information" would redesign an excitation signal while the real fault is a
/// covariance that stopped being a covariance three samples ago.
///
///  * non_finite_state            : the carried covariance or parameter vector is
///                                  already non-finite when the cycle begins.
///                                  Nothing this cycle can produce is
///                                  meaningful, and the recursion has no
///                                  mechanism that returns a non-finite
///                                  covariance to a finite one, so only
///                                  reconstruction recovers.
///  * non_finite_observation      : the supplied output sample is not finite. It
///                                  enters the prediction error and from there
///                                  the parameter step directly, so admitting
///                                  one destroys the parameters permanently
///                                  while every downstream query still reports
///                                  numbers. The repair is upstream, in whatever
///                                  measures the output.
///  * non_finite_regressor        : the supplied regressor has a non-finite
///                                  component. It reaches the parameters through
///                                  the gain and the covariance through the rank
///                                  one update, so it destroys both. A different
///                                  subsystem from the one that measures the
///                                  output, hence a different enumerator.
///  * non_finite_denominator      : every operand was finite and the
///                                  covariance-weighted regressor still left the
///                                  scalar's range. That is a magnitude fault in
///                                  the carried covariance or the regressor
///                                  scale, not a domain violation of either, and
///                                  it is worth its own name because the repair
///                                  is a rescaling rather than a replacement.
///  * indefinite_covariance       : the denominator is negative by more than the
///                                  resolution below. It is lambda + phi' P phi
///                                  with lambda strictly positive, so a negative
///                                  value means P is no longer positive
///                                  semidefinite and the gain formed from it
///                                  would point AGAINST the prediction error --
///                                  the same failure the negative forgetting
///                                  factor is rejected for, arrived at through
///                                  the recursion instead of through the
///                                  configuration.
///  * denominator_below_resolution: the denominator carries no significant
///                                  digits at the scale of the operands that
///                                  formed it. Two situations reach this and it
///                                  does NOT distinguish them, deliberately,
///                                  because the arithmetic cannot: a regressor
///                                  the covariance genuinely cannot resolve, and
///                                  a cancellation in phi' P phi that leaves a
///                                  residue made entirely of rounding. A gain
///                                  divided by either is noise.
enum class rls_update_error
{
    non_finite_state,
    non_finite_observation,
    non_finite_regressor,
    non_finite_denominator,
    indefinite_covariance,
    denominator_below_resolution,
    non_finite_result,
};

template <typename Scalar, std::size_t NP>
struct rls_config
{
    Scalar lambda{Scalar{0.99}};
    Matrix<Scalar, NP, NP> P0{Matrix<Scalar, NP, NP>::Identity() * Scalar{1000}};
    Scalar cov_upper_bound{Scalar{1e6}};
};

template <ctrlpp_floating_scalar Scalar, std::size_t NP>
class rls
{
    static_assert(NP >= 1, "rls requires at least one parameter");

public:
    /// @brief Fallible factory, and the only way to originate an estimator.
    ///
    /// The forgetting factor's condition is derived from the covariance update
    /// below rather than asserted as a range: lambda is the divisor of that
    /// update, which fixes finite-and-positive as the exact domain. The upper
    /// half of the documented (0, 1] contract is enforced separately and named
    /// as a contract violation, because the arithmetic there is well defined and
    /// its consequence -- a covariance deflated until the estimator stops
    /// adapting -- is silent rather than explosive. See `rls_error`.
    static auto create(rls_config<Scalar, NP> config = {}) -> ctrlpp::expected<rls, rls_error>
    {
        if(!std::isfinite(config.lambda) || !(config.lambda > Scalar{0}))
            return ctrlpp::unexpected(rls_error::non_positive_forgetting_factor);
        if(config.lambda > Scalar{1})
            return ctrlpp::unexpected(rls_error::forgetting_factor_above_unity);
        if(!config.P0.allFinite())
            return ctrlpp::unexpected(rls_error::non_finite_initial_covariance);
        if(!std::isfinite(config.cov_upper_bound) || !(config.cov_upper_bound > Scalar{0}))
            return ctrlpp::unexpected(rls_error::non_positive_covariance_bound);
        return rls{validated_tag{}, std::move(config)};
    }

    /// @brief Rounded operations behind the update's denominator.
    ///
    /// Enumerated rather than chosen. The denominator is
    /// lambda + phi' (P phi), and two contractions over the parameter dimension
    /// stand behind it: forming P phi costs NP multiplies and NP - 1 additions
    /// per component, and contracting phi against it costs the same again, that
    /// is 2 * (2 * NP - 1); the final sum with the forgetting factor is one more.
    /// Every operation is counted whether or not it rounds, so the count bounds
    /// the accumulated error from above rather than describing it tightly, which
    /// is what a resolution floor requires.
    static constexpr int denominator_rounding_ops = 2 * (2 * static_cast<int>(NP) - 1) + 1;

    /// @brief Incorporate one observation, or report why the cycle was refused.
    ///
    /// The cycle is classified before any member is written, so a refused cycle
    /// leaves the parameters and the covariance bitwise unchanged and the caller
    /// may retry on the next sample. That matters more here than the boolean it
    /// replaces suggested: the parameters and the covariance are the estimator's
    /// entire memory, nothing re-derives them, so one admitted non-finite sample
    /// makes them non-finite forever while `parameters()` keeps returning a
    /// vector the caller has no way to distrust.
    ///
    /// **The denominator's floor is derived, not chosen.** It is the counted
    /// rounding of the two contractions that formed it, at the scale of the
    /// largest operand that entered -- the forgetting factor, or the
    /// Cauchy-Schwarz bound on the quadratic form. An absolute floor is wrong in
    /// both directions and both are reachable: on a problem posed far below unit
    /// scale it refuses a perfectly well conditioned update, and on one posed far
    /// above it accepts a denominator whose significant digits have all canceled
    /// away and divides a gain by rounding noise.
    auto update(Scalar y, const Vector<Scalar, NP>& phi) -> ctrlpp::expected<void, rls_update_error>
    {
        if(!m_P.allFinite() || !m_theta.allFinite())
            return ctrlpp::unexpected(rls_update_error::non_finite_state);
        if(!std::isfinite(y))
            return ctrlpp::unexpected(rls_update_error::non_finite_observation);
        if(!phi.allFinite())
            return ctrlpp::unexpected(rls_update_error::non_finite_regressor);

        Vector<Scalar, NP> P_phi = m_P * phi;
        Scalar denom = m_lambda + phi.dot(P_phi);

        if(!std::isfinite(denom))
            return ctrlpp::unexpected(rls_update_error::non_finite_denominator);

        Scalar const relative_floor = static_cast<Scalar>(denominator_rounding_ops)
                                      * std::numeric_limits<Scalar>::epsilon();
        Scalar const abs_denom = std::abs(denom);
        Scalar const phi_norm = phi.stableNorm();
        Scalar const P_phi_norm = P_phi.stableNorm();

        if(at_or_below_product(abs_denom, relative_floor, m_lambda, Scalar{1})
           || at_or_below_product(abs_denom, relative_floor, phi_norm, P_phi_norm))
            return ctrlpp::unexpected(rls_update_error::denominator_below_resolution);
        if(denom < Scalar{0})
            return ctrlpp::unexpected(rls_update_error::indefinite_covariance);

        Scalar e = y - phi.dot(m_theta);
        Vector<Scalar, NP> k = P_phi / denom;

        auto next_theta = (m_theta + k * e).eval();
        auto next_P =
            ((m_P - k * P_phi.transpose()).eval() / m_lambda).eval();
        next_P = ((next_P + next_P.transpose()) * Scalar{0.5}).eval();

        Scalar trace = next_P.trace();
        Scalar trace_bound = m_cov_upper_bound * static_cast<Scalar>(NP);
        if(trace > trace_bound)
            next_P *= trace_bound / trace;

        if(!std::isfinite(e) || !k.allFinite() || !next_theta.allFinite()
            || !next_P.allFinite() || !std::isfinite(trace)
            || !std::isfinite(trace_bound))
            return ctrlpp::unexpected(rls_update_error::non_finite_result);

        m_theta = std::move(next_theta);
        m_P = std::move(next_P);
        return {};
    }

    const Vector<Scalar, NP>& parameters() const { return m_theta; }

    const Matrix<Scalar, NP, NP>& covariance() const { return m_P; }

private:
    static auto at_or_below_product(Scalar value, Scalar factor, Scalar lhs, Scalar rhs) -> bool
    {
        if(lhs == Scalar{0} || rhs == Scalar{0})
            return value == Scalar{0};

        int value_exponent{};
        int factor_exponent{};
        int lhs_exponent{};
        int rhs_exponent{};
        int product_exponent{};

        Scalar const value_fraction = std::frexp(value, &value_exponent);
        Scalar const factor_fraction = std::frexp(factor, &factor_exponent);
        Scalar const lhs_fraction = std::frexp(lhs, &lhs_exponent);
        Scalar const rhs_fraction = std::frexp(rhs, &rhs_exponent);
        Scalar const product_fraction =
            std::frexp(factor_fraction * lhs_fraction * rhs_fraction, &product_exponent);

        int const threshold_exponent =
            factor_exponent + lhs_exponent + rhs_exponent + product_exponent;
        if(value_exponent != threshold_exponent)
            return value_exponent < threshold_exponent;
        return value_fraction <= product_fraction;
    }

    /// @brief Tag selecting the non-validating constructor reserved for
    /// `create`, which is what makes the factory the only public path and the
    /// validation impossible to bypass.
    struct validated_tag
    {
    };

    rls(validated_tag, rls_config<Scalar, NP> config) : m_lambda{config.lambda}, m_cov_upper_bound{config.cov_upper_bound}, m_theta{Vector<Scalar, NP>::Zero()}, m_P{std::move(config.P0)} {}

    Scalar m_lambda;
    Scalar m_cov_upper_bound;
    Vector<Scalar, NP> m_theta;
    Matrix<Scalar, NP, NP> m_P;
};

}

#endif
