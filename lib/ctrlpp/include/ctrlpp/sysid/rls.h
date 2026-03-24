#ifndef HPP_GUARD_CTRLPP_SYSID_RLS_H
#define HPP_GUARD_CTRLPP_SYSID_RLS_H

/// @brief Recursive Least Squares with bounded covariance and forgetting factor.
///
/// @cite ljung1999 -- Ljung, "System Identification: Theory for the User", 1999

#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <cstddef>

namespace ctrlpp
{

template <typename Scalar, std::size_t NP>
struct rls_config
{
    Scalar lambda{Scalar{0.99}};
    Matrix<Scalar, NP, NP> P0{Matrix<Scalar, NP, NP>::Identity() * Scalar{1000}};
    Scalar cov_upper_bound{Scalar{1e6}};
};

template <typename Scalar, std::size_t NP>
class rls
{
public:
    explicit rls(rls_config<Scalar, NP> config = {}) : m_lambda{config.lambda}, m_cov_upper_bound{config.cov_upper_bound}, m_theta{Vector<Scalar, NP>::Zero()}, m_P{config.P0} {}

    void update(Scalar y, const Vector<Scalar, NP>& phi)
    {
        Scalar e = y - phi.dot(m_theta);

        Vector<Scalar, NP> P_phi = m_P * phi;
        Scalar denom = m_lambda + phi.dot(P_phi);
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
    }

    const Vector<Scalar, NP>& parameters() const { return m_theta; }

    const Matrix<Scalar, NP, NP>& covariance() const { return m_P; }

private:
    Scalar m_lambda;
    Scalar m_cov_upper_bound;
    Vector<Scalar, NP> m_theta;
    Matrix<Scalar, NP, NP> m_P;
};

} // namespace ctrlpp

#endif
