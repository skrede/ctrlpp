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
#include <cstddef>
#include <utility>
#include <algorithm>

namespace ctrlpp
{

template <ctrlpp_floating_scalar Scalar, std::size_t NA, std::size_t NB, std::size_t NU = 1, std::size_t NY = 1>
class recursive_arx
{
    static_assert(NA >= 1 && NB >= 1, "recursive_arx requires NA >= 1 and NB >= 1");

public:
    static constexpr std::size_t NP = NA * NY + NB * NU;

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
    /// The refusal is the estimator's own, forwarded verbatim rather than
    /// restated under a second name that could drift out of step with the
    /// arithmetic it describes. See `rls_update_error`.
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
    auto update(Scalar y, Scalar u) -> ctrlpp::expected<void, rls_update_error>
    {
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
            return ctrlpp::unexpected(applied.error());

        m_y_hist[m_write_idx % NA] = y;
        m_u_hist[m_write_idx % NB] = u;
        ++m_write_idx;
        ++m_sample_count;
        return {};
    }

    const Vector<Scalar, NP>& parameters() const { return m_rls.parameters(); }

    const Matrix<Scalar, NP, NP>& covariance() const { return m_rls.covariance(); }

    discrete_state_space<Scalar, std::max(NA, NB), NU, NY> to_state_space() const
    {
        auto theta = m_rls.parameters();

        // Observer canonical realization dimension = max(deg A, deg B) (Ljung 1999, Ch. 4).
        // When NB > NA the extra b-coefficients b_{NA+1..NB} require additional states.
        static constexpr std::size_t NX = std::max(NA, NB);

        Matrix<Scalar, NX, NX> A = Matrix<Scalar, NX, NX>::Zero();
        Matrix<Scalar, NX, NU> B = Matrix<Scalar, NX, NU>::Zero();
        Matrix<Scalar, NY, NX> C = Matrix<Scalar, NY, NX>::Zero();
        Matrix<Scalar, NY, NU> D = Matrix<Scalar, NY, NU>::Zero();

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
