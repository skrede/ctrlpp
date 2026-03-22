#ifndef HPP_GUARD_CTRLPP_SYSID_RLS_H
#define HPP_GUARD_CTRLPP_SYSID_RLS_H

#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <cstddef>

namespace ctrlpp {

template<typename Scalar, std::size_t NP>
struct rls_config {
    Scalar lambda{Scalar{0.99}};
    Matrix<Scalar, NP, NP> P0{Matrix<Scalar, NP, NP>::Identity() * Scalar{1000}};
    Scalar cov_upper_bound{Scalar{1e6}};
};

template<typename Scalar, std::size_t NP>
class rls {
  public:
    explicit rls(rls_config<Scalar, NP> config = {})
        : lambda_{config.lambda}
        , cov_upper_bound_{config.cov_upper_bound}
        , theta_{Vector<Scalar, NP>::Zero()}
        , P_{config.P0}
    {
    }

    void update(Scalar y, const Vector<Scalar, NP>& phi)
    {
        Scalar e = y - phi.dot(theta_);

        Vector<Scalar, NP> P_phi = P_ * phi;
        Scalar denom = lambda_ + phi.dot(P_phi);
        Vector<Scalar, NP> k = P_phi / denom;

        theta_ += k * e;

        P_ = (P_ - k * P_phi.transpose()).eval() / lambda_;
        P_ = (P_ + P_.transpose()) * Scalar{0.5};

        Scalar trace = P_.trace();
        Scalar trace_bound = cov_upper_bound_ * static_cast<Scalar>(NP);
        if (trace > trace_bound) {
            P_ *= trace_bound / trace;
        }
    }

    [[nodiscard]] auto parameters() const -> const Vector<Scalar, NP>&
    {
        return theta_;
    }

    [[nodiscard]] auto covariance() const -> const Matrix<Scalar, NP, NP>&
    {
        return P_;
    }

  private:
    Scalar lambda_;
    Scalar cov_upper_bound_;
    Vector<Scalar, NP> theta_;
    Matrix<Scalar, NP, NP> P_;
};

}

#endif
