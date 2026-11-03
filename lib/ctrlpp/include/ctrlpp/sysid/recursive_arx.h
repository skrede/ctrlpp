#ifndef HPP_GUARD_CTRLPP_SYSID_RECURSIVE_ARX_H
#define HPP_GUARD_CTRLPP_SYSID_RECURSIVE_ARX_H

// @cite ljung1999 -- Ljung, "System Identification: Theory for the User", 1999, Ch. 11

#include "ctrlpp/sysid/rls.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

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

    // Only create() can name the tag, so the estimator's validation cannot be bypassed.
    struct validated_tag
    {
    };

public:
    static constexpr std::size_t NP = NA + NB;

    static ctrlpp::expected<recursive_arx, rls_error> create(rls_config<Scalar, NP> config = {})
    {
        auto estimator = rls<Scalar, NP>::create(std::move(config));
        if(!estimator)
            return ctrlpp::unexpected(estimator.error());
        return recursive_arx{validated_tag{}, std::move(*estimator)};
    }

    // The history is written only once the estimator accepts the sample, since it is the next cycle's regressor.
    ctrlpp::expected<void, recursive_arx_update_error> update(Scalar y, Scalar u)
    {
        if(!std::isfinite(u))
            return ctrlpp::unexpected(recursive_arx_update_error::non_finite_input);

        Vector<Scalar, NP> phi = Vector<Scalar, NP>::Zero();

        for(std::size_t i = 0; i < NA; ++i)
        {
            std::size_t idx = (m_write_idx + NA - 1 - i) % NA;
            phi(static_cast<Eigen::Index>(i)) = m_y_hist[idx];
        }
        for(std::size_t i = 0; i < NB; ++i)
        {
            std::size_t idx = (m_write_idx + NB - 1 - i) % NB;
            phi(static_cast<Eigen::Index>(NA + i)) = m_u_hist[idx];
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

        // Observer canonical realization of dimension max(deg A, deg B) (Ljung 1999, Ch. 4): when NB > NA the extra
        // b-coefficients need states of their own.
        static constexpr std::size_t NX = std::max(NA, NB);

        Matrix<Scalar, NX, NX> A = Matrix<Scalar, NX, NX>::Zero();
        Matrix<Scalar, NX, 1> B = Matrix<Scalar, NX, 1>::Zero();
        Matrix<Scalar, 1, NX> C = Matrix<Scalar, 1, NX>::Zero();
        Matrix<Scalar, 1, 1> D = Matrix<Scalar, 1, 1>::Zero();

        for(std::size_t i = 0; i < NA; ++i)
            A(static_cast<Eigen::Index>(i), 0) = theta(static_cast<Eigen::Index>(i));
        for(std::size_t i = 0; i + 1 < NX; ++i)
            A(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(i + 1)) = Scalar{1};

        for(std::size_t i = 0; i < NB; ++i)
            B(static_cast<Eigen::Index>(i), 0) = theta(static_cast<Eigen::Index>(NA + i));

        C(0, 0) = Scalar{1};

        return {A, B, C, D};
    }

private:
    rls<Scalar, NP> m_rls;
    std::size_t m_write_idx;
    std::size_t m_sample_count;
    std::array<Scalar, NA> m_y_hist;
    std::array<Scalar, NB> m_u_hist;

    recursive_arx(validated_tag, rls<Scalar, NP> estimator)
        : m_rls{std::move(estimator)}
        , m_write_idx{0}
        , m_sample_count{0}
        , m_y_hist{}
        , m_u_hist{}
    {
    }

    static constexpr recursive_arx_update_error forward_error(rls_update_error error)
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
};

}

#endif
