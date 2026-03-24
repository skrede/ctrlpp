#ifndef HPP_GUARD_CTRLPP_SYSID_RECURSIVE_ARX_H
#define HPP_GUARD_CTRLPP_SYSID_RECURSIVE_ARX_H

#include "ctrlpp/types.h"
#include "ctrlpp/model/state_space.h"

#include "ctrlpp/sysid/rls.h"

#include <Eigen/Dense>

#include <array>
#include <cstddef>

namespace ctrlpp {

template <typename Scalar, std::size_t NA, std::size_t NB, std::size_t NU = 1, std::size_t NY = 1>
class recursive_arx
{
public:
    static constexpr std::size_t NP = NA * NY + NB * NU;

    explicit recursive_arx(rls_config<Scalar, NP> config = {})
        : m_rls{config}
    {
    }

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

    const Vector<Scalar, NP> &parameters() const
    {
        return m_rls.parameters();
    }

    const Matrix<Scalar, NP, NP> &covariance() const
    {
        return m_rls.covariance();
    }

    discrete_state_space<Scalar, NA, NU, NY> to_state_space() const
    {
        auto theta = m_rls.parameters();

        Matrix<Scalar, NA, NA> A = Matrix<Scalar, NA, NA>::Zero();
        Matrix<Scalar, NA, NU> B = Matrix<Scalar, NA, NU>::Zero();
        Matrix<Scalar, NY, NA> C = Matrix<Scalar, NY, NA>::Zero();
        Matrix<Scalar, NY, NU> D = Matrix<Scalar, NY, NU>::Zero();

        // Observer canonical form for ARX:
        //   y(t) = a1*y(t-1) + ... + aNa*y(t-NA) + b1*u(t-1) + ... + bNb*u(t-NB)
        //
        //   A = [a1  1  0 ...]    B = [b1]    C = [1 0 ... 0]    D = [0]
        //       [a2  0  1 ...]        [b2]
        //       [... ... ...]         [...]
        //       [aNa 0  0 ...]        [bNa (or 0)]
        //
        // A: first column = a-coefficients, superdiagonal = 1
        for(std::size_t i = 0; i < NA; ++i)
            A(static_cast<int>(i), 0) = theta(static_cast<int>(i));
        for(std::size_t i = 0; i + 1 < NA; ++i)
            A(static_cast<int>(i), static_cast<int>(i + 1)) = Scalar{1};

        // B column: b-coefficients, zero-padded if NB < NA
        for(std::size_t i = 0; i < NB && i < NA; ++i)
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
