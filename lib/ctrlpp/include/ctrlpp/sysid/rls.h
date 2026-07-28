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
#include <cstddef>
#include <utility>

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

    /// @brief Incorporate one observation.
    ///
    /// Returns whether the update was applied. That boolean is a narrower report
    /// than the failure channel this type's construction now uses -- it cannot
    /// say why a sample was skipped -- and converting it is separately owned
    /// work; the contract is unchanged here.
    bool update(Scalar y, const Vector<Scalar, NP>& phi)
    {
        Scalar e = y - phi.dot(m_theta);

        Vector<Scalar, NP> P_phi = m_P * phi;
        Scalar denom = m_lambda + phi.dot(P_phi);
        if(!std::isfinite(denom) || std::abs(denom) < Scalar{1e-14})
            return false;
        Vector<Scalar, NP> k = P_phi / denom;

        m_theta += k * e;

        m_P = (m_P - k * P_phi.transpose()).eval() / m_lambda;
        m_P = (m_P + m_P.transpose()) * Scalar{0.5};

        Scalar trace = m_P.trace();
        Scalar trace_bound = m_cov_upper_bound * static_cast<Scalar>(NP);
        if(trace > trace_bound)
        {
            m_P *= trace_bound / trace;
        }
        return true;
    }

    const Vector<Scalar, NP>& parameters() const { return m_theta; }

    const Matrix<Scalar, NP, NP>& covariance() const { return m_P; }

private:
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
