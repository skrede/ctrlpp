#ifndef HPP_GUARD_CTRLPP_SYSID_RLS_H
#define HPP_GUARD_CTRLPP_SYSID_RLS_H

// @cite ljung1999 -- Ljung, "System Identification: Theory for the User", 1999

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <algorithm>

namespace ctrlpp
{

enum class rls_error
{
    non_positive_forgetting_factor,
    forgetting_factor_above_unity,
    non_finite_initial_covariance,
    non_positive_covariance_bound,
};

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

    // Only create() can name the tag, so its validation cannot be bypassed.
    struct validated_tag
    {
    };

public:
    // Counted, not chosen: P phi and phi' (P phi) each cost 2 * NP - 1 operations, and adding lambda one more.
    static constexpr std::int32_t denominator_rounding_ops = 2 * (2 * static_cast<std::int32_t>(NP) - 1) + 1;

    static ctrlpp::expected<rls, rls_error> create(rls_config<Scalar, NP> config = {})
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

    // Every refusal is decided before a member is written, so a refused cycle leaves the estimator bitwise unchanged.
    ctrlpp::expected<void, rls_update_error> update(Scalar y, const Vector<Scalar, NP>& phi)
    {
        if(!m_P.allFinite() || !m_theta.allFinite())
            return ctrlpp::unexpected(rls_update_error::non_finite_state);
        if(!std::isfinite(y))
            return ctrlpp::unexpected(rls_update_error::non_finite_observation);
        if(!phi.allFinite())
            return ctrlpp::unexpected(rls_update_error::non_finite_regressor);

        Vector<Scalar, NP> P_phi = m_P * phi;
        Scalar denom = m_lambda + phi.dot(P_phi);
        if(auto const resolved = check_denominator(denom, phi, P_phi); !resolved)
            return resolved;
        return apply_update(y, phi, P_phi, denom);
    }

    const Vector<Scalar, NP>& parameters() const { return m_theta; }

    const Matrix<Scalar, NP, NP>& covariance() const { return m_P; }

private:
    Scalar m_lambda;
    Scalar m_cov_upper_bound;
    Vector<Scalar, NP> m_theta;
    Matrix<Scalar, NP, NP> m_P;

    rls(validated_tag, rls_config<Scalar, NP> config)
        : m_lambda{config.lambda}
        , m_cov_upper_bound{config.cov_upper_bound}
        , m_theta{Vector<Scalar, NP>::Zero()}
        , m_P{std::move(config.P0)}
    {
    }

    ctrlpp::expected<void, rls_update_error> check_denominator(Scalar denom, const Vector<Scalar, NP>& phi, const Vector<Scalar, NP>& P_phi) const
    {
        if(!std::isfinite(denom))
            return ctrlpp::unexpected(rls_update_error::non_finite_denominator);

        Scalar const relative_floor = static_cast<Scalar>(denominator_rounding_ops) * std::numeric_limits<Scalar>::epsilon();
        Scalar const abs_denom = std::abs(denom);
        Scalar const phi_norm = phi.stableNorm();
        Scalar const P_phi_norm = P_phi.stableNorm();

        if(at_or_below_product(abs_denom, relative_floor, m_lambda, Scalar{1}) || at_or_below_product(abs_denom, relative_floor, phi_norm, P_phi_norm))
            return ctrlpp::unexpected(rls_update_error::denominator_below_resolution);
        if(denom < Scalar{0})
            return ctrlpp::unexpected(rls_update_error::indefinite_covariance);
        return {};
    }

    ctrlpp::expected<void, rls_update_error> apply_update(Scalar y, const Vector<Scalar, NP>& phi, const Vector<Scalar, NP>& P_phi, Scalar denom)
    {
        Scalar e = y - phi.dot(m_theta);
        Vector<Scalar, NP> k = P_phi / denom;

        auto next_theta = (m_theta + k * e).eval();
        auto next_P = ((m_P - k * P_phi.transpose()).eval() / m_lambda).eval();
        next_P = ((next_P + next_P.transpose()) * Scalar{0.5}).eval();

        Scalar trace = next_P.trace();
        Scalar trace_bound = m_cov_upper_bound * static_cast<Scalar>(NP);
        if(trace > trace_bound)
            next_P *= trace_bound / trace;

        if(!std::isfinite(e) || !k.allFinite() || !next_theta.allFinite() || !next_P.allFinite() || !std::isfinite(trace) || !std::isfinite(trace_bound))
            return ctrlpp::unexpected(rls_update_error::non_finite_result);

        m_theta = std::move(next_theta);
        m_P = std::move(next_P);
        return {};
    }

    static bool at_or_below_product(Scalar value, Scalar factor, Scalar lhs, Scalar rhs)
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
        Scalar const product_fraction = std::frexp(factor_fraction * lhs_fraction * rhs_fraction, &product_exponent);

        int const threshold_exponent = factor_exponent + lhs_exponent + rhs_exponent + product_exponent;
        if(value_exponent != threshold_exponent)
            return value_exponent < threshold_exponent;
        return value_fraction <= product_fraction;
    }
};

}

#endif
