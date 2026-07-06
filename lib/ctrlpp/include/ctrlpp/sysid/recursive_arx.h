#ifndef HPP_GUARD_CTRLPP_SYSID_RECURSIVE_ARX_H
#define HPP_GUARD_CTRLPP_SYSID_RECURSIVE_ARX_H

/// @brief Recursive ARX model identification using RLS with state-space conversion.
///
/// @cite ljung1999 -- Ljung, "System Identification: Theory for the User", 1999, Ch. 11

#include "ctrlpp/types.h"

#include "ctrlpp/sysid/rls.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/model/state_space.h"

#include <Eigen/Dense>

#include <array>
#include <cstddef>
#include <algorithm>

namespace ctrlpp
{

template <ctrlpp_floating_scalar Scalar, std::size_t NA, std::size_t NB, std::size_t NU = 1, std::size_t NY = 1>
class recursive_arx
{
    static_assert(NA >= 1 && NB >= 1, "recursive_arx requires NA >= 1 and NB >= 1");

public:
    static constexpr std::size_t NP = NA * NY + NB * NU;

    explicit recursive_arx(rls_config<Scalar, NP> config = {}) : m_rls{config} {}

    void update(Scalar y, Scalar u)
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

        m_rls.update(y, phi);

        m_y_hist[m_write_idx % NA] = y;
        m_u_hist[m_write_idx % NB] = u;
        ++m_write_idx;
        ++m_sample_count;
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
    rls<Scalar, NP> m_rls;
    std::array<Scalar, NA> m_y_hist{};
    std::array<Scalar, NB> m_u_hist{};
    std::size_t m_write_idx{0};
    std::size_t m_sample_count{0};
};

}

#endif
