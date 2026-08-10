#ifndef HPP_GUARD_CTRLPP_SYSID_RECURSIVE_ARX_H
#define HPP_GUARD_CTRLPP_SYSID_RECURSIVE_ARX_H

/// @brief Recursive ARX model identification using RLS with state-space conversion.
///
/// @cite ljung1999 -- Ljung, "System Identification: Theory for the User", 1999, Ch. 11

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/sysid/rls.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/model/state_space.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <cstddef>
#include <utility>
#include <algorithm>

namespace ctrlpp
{

enum class recursive_arx_update_error
{
    non_finite_input,
    non_finite_state,
    non_finite_observation,
    non_finite_regressor,
    non_finite_denominator,
    indefinite_covariance,
    denominator_below_resolution,
    non_finite_result,
};

template <ctrlpp_floating_scalar Scalar, std::size_t NA, std::size_t NB>
class recursive_arx
{
    static_assert(NA >= 1 && NB >= 1, "recursive_arx requires NA >= 1 and NB >= 1");

public:
    static constexpr std::size_t NP = NA + NB;

    /// @brief Fallible factory, and the only way to originate an estimator.
    ///
    /// This type owns a recursive least-squares estimator and configures it from
    /// the same aggregate, so it has no configuration condition of its own to
    /// state: it forwards that estimator's rejection verbatim rather than
    /// restating the conditions here, where a second copy could drift out of
    /// step with the arithmetic it describes. See `rls_error`.
    static auto create(rls_config<Scalar, NP> config = {}) -> ctrlpp::expected<recursive_arx, rls_error>
    {
        auto estimator = rls<Scalar, NP>::create(std::move(config));
        if(!estimator)
            return ctrlpp::unexpected(estimator.error());
        return recursive_arx{validated_tag{}, std::move(*estimator)};
    }

    /// @brief Incorporate one input-output sample, or report why the cycle was
    /// refused.
    ///
    /// This wrapper previously called the estimator and discarded its answer, so
    /// a refused sample was swallowed here and no caller could learn that the
    /// model had stopped moving.
    ///
    /// A refused cycle also leaves the regressor history, the write index and the
    /// sample count untouched. That is not tidiness: the history IS the next
    /// cycle's regressor, so recording a sample the estimator refused as
    /// non-finite would poison every regressor built afterwards -- the poison
    /// would latch in the wrapper after the estimator had correctly declined it.
    auto update(Scalar y, Scalar u) -> ctrlpp::expected<void, recursive_arx_update_error>
    {
        if(!std::isfinite(u))
            return ctrlpp::unexpected(recursive_arx_update_error::non_finite_input);

        Vector<Scalar, NP> phi = Vector<Scalar, NP>::Zero();

        // Build regressor: [y(t-1), ..., y(t-NA), u(t-1), ..., u(t-NB)]
        for(std::size_t i = 0; i < NA; ++i)
        {
            std::size_t idx = (m_write_idx + NA - 1 - i) % NA;
            phi(static_cast<int>(i)) = m_y_hist[idx];
        }
        for(std::size_t i = 0; i < NB; ++i)
        {
            std::size_t idx = (m_write_idx + NB - 1 - i) % NB;
            phi(static_cast<int>(NA + i)) = m_u_hist[idx];
        }

        if(const auto applied = m_rls.update(y, phi); !applied)
            return ctrlpp::unexpected(forward_error(applied.error()));

        m_y_hist[m_write_idx % NA] = y;
        m_u_hist[m_write_idx % NB] = u;
        ++m_write_idx;
        ++m_sample_count;
        return {};
    }

    const Vector<Scalar, NP>& parameters() const { return m_rls.parameters(); }

    const Matrix<Scalar, NP, NP>& covariance() const { return m_rls.covariance(); }

    discrete_state_space<Scalar, std::max(NA, NB), 1, 1> to_state_space() const
    {
        auto theta = m_rls.parameters();

        // Observer canonical realization dimension = max(deg A, deg B) (Ljung 1999, Ch. 4).
        // When NB > NA the extra b-coefficients b_{NA+1..NB} require additional states.
        static constexpr std::size_t NX = std::max(NA, NB);

        Matrix<Scalar, NX, NX> A = Matrix<Scalar, NX, NX>::Zero();
        Matrix<Scalar, NX, 1> B = Matrix<Scalar, NX, 1>::Zero();
        Matrix<Scalar, 1, NX> C = Matrix<Scalar, 1, NX>::Zero();
        Matrix<Scalar, 1, 1> D = Matrix<Scalar, 1, 1>::Zero();

        // Observer canonical form for ARX (NX = max(NA, NB) states):
        //   y(t) = a1*y(t-1) + ... + aNa*y(t-NA) + b1*u(t-1) + ... + bNb*u(t-NB)
        //
        //   A = [a1  1  0 ...]    B = [b1]    C = [1 0 ... 0]    D = [0]
        //       [a2  0  1 ...]        [b2]
        //       [... ... ...]         [...]
        //       [ .  0  0 ...]        [bNx (or 0)]
        //
        // A: first column = a-coefficients (rows NA..NX-1 stay zero), superdiagonal = 1
        for(std::size_t i = 0; i < NA; ++i)
            A(static_cast<int>(i), 0) = theta(static_cast<int>(i));
        for(std::size_t i = 0; i + 1 < NX; ++i)
            A(static_cast<int>(i), static_cast<int>(i + 1)) = Scalar{1};

        // B column: all b-coefficients (rows NB..NX-1 stay zero)
        for(std::size_t i = 0; i < NB; ++i)
            B(static_cast<int>(i), 0) = theta(static_cast<int>(NA + i));

        // C = [1, 0, ..., 0]
        C(0, 0) = Scalar{1};

        return {.A = A, .B = B, .C = C, .D = D};
    }

private:
    static constexpr auto forward_error(rls_update_error error) -> recursive_arx_update_error
    {
        switch(error)
        {
            case rls_update_error::non_finite_state:
                return recursive_arx_update_error::non_finite_state;
            case rls_update_error::non_finite_observation:
                return recursive_arx_update_error::non_finite_observation;
            case rls_update_error::non_finite_regressor:
                return recursive_arx_update_error::non_finite_regressor;
            case rls_update_error::non_finite_denominator:
                return recursive_arx_update_error::non_finite_denominator;
            case rls_update_error::indefinite_covariance:
                return recursive_arx_update_error::indefinite_covariance;
            case rls_update_error::denominator_below_resolution:
                return recursive_arx_update_error::denominator_below_resolution;
            case rls_update_error::non_finite_result:
                return recursive_arx_update_error::non_finite_result;
        }
        return recursive_arx_update_error::non_finite_state;
    }

    /// @brief Tag selecting the non-validating constructor reserved for
    /// `create`, which is what makes the factory the only public path and the
    /// validation impossible to bypass.
    struct validated_tag
    {
    };

    recursive_arx(validated_tag, rls<Scalar, NP> estimator) : m_rls{std::move(estimator)} {}

    rls<Scalar, NP> m_rls;
    std::array<Scalar, NA> m_y_hist{};
    std::array<Scalar, NB> m_u_hist{};
    std::size_t m_write_idx{0};
    std::size_t m_sample_count{0};
};

}

#endif
