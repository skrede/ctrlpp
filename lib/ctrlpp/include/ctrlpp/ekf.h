#ifndef HPP_GUARD_CTRLPP_EKF_H
#define HPP_GUARD_CTRLPP_EKF_H

#include "ctrlpp/types.h"
#include "ctrlpp/observer_policy.h"

#include "ctrlpp/detail/numerical_diff.h"

#include "ctrlpp/mpc/dynamics_model.h"
#include "ctrlpp/mpc/measurement_model.h"
#include "ctrlpp/mpc/differentiable_dynamics.h"
#include "ctrlpp/mpc/differentiable_measurement.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <utility>

namespace ctrlpp {

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct ekf_config
{
    Matrix<Scalar, NX, NX> Q{Matrix<Scalar, NX, NX>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Vector<Scalar, NX> x0{Vector<Scalar, NX>::Zero()};
    Matrix<Scalar, NX, NX> P0{Matrix<Scalar, NX, NX>::Identity()};
    Scalar numerical_eps{std::sqrt(std::numeric_limits<Scalar>::epsilon())};
};

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, typename Dynamics, typename Measurement>
    requires dynamics_model<Dynamics, Scalar, NX, NU> && measurement_model<Measurement, Scalar, NX, NY>
class ekf
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);
    static constexpr int ny = static_cast<int>(NY);

public:
    using observer_tag = struct ekf_tag;
    using state_vector_t = Vector<Scalar, NX>;
    using input_vector_t = Vector<Scalar, NU>;
    using output_vector_t = Vector<Scalar, NY>;
    using cov_matrix_t = Matrix<Scalar, NX, NX>;
    using meas_cov_matrix_t = Matrix<Scalar, NY, NY>;

    ekf(Dynamics dynamics, Measurement measurement, ekf_config<Scalar, NX, NU, NY> config)
        : m_eps{config.numerical_eps}
        , m_dynamics{std::move(dynamics)}
        , m_P{std::move(config.P0)}
        , m_Q{std::move(config.Q)}
        , m_x{std::move(config.x0)}
        , m_R{std::move(config.R)}
        , m_measurement{std::move(measurement)}
        , m_innovation{output_vector_t::Zero()}
    {
    }

    void predict(const input_vector_t &u)
    {
        auto x_prev = m_x;

        m_x = m_dynamics(x_prev, u);

        Matrix<Scalar, NX, NX> F;
        if constexpr(differentiable_dynamics<Dynamics, Scalar, NX, NU>)
            F = m_dynamics.jacobian_x(x_prev, u);
        else
            F = detail::numerical_jacobian_x<Scalar, NX, NU>(m_dynamics, x_prev, u, m_eps);

        m_P = (F * m_P * F.transpose() + m_Q).eval();
        m_P = (Scalar{0.5} * (m_P + m_P.transpose())).eval();
    }

    void update(const output_vector_t &z)
    {
        auto z_pred = m_measurement(m_x);

        m_innovation = (z - z_pred).eval();

        Matrix<Scalar, NY, NX> H;
        if constexpr(differentiable_measurement<Measurement, Scalar, NX, NY>)
            H = m_measurement.jacobian(m_x);
        else
            H = detail::numerical_jacobian_h<Scalar, NX, NY>(m_measurement, m_x, m_eps);

        meas_cov_matrix_t S = (H * m_P * H.transpose() + m_R).eval();

        Eigen::Matrix<Scalar, ny, nx> KT_solved = S.transpose().colPivHouseholderQr().solve(H * m_P);
        Eigen::Matrix<Scalar, nx, ny> K = KT_solved.transpose().eval();

        m_x = (m_x + K * m_innovation).eval();

        cov_matrix_t IKH = cov_matrix_t::Identity() - K * H;
        m_P = (IKH * m_P * IKH.transpose() + K * m_R * K.transpose()).eval();
        m_P = (Scalar{0.5} * (m_P + m_P.transpose())).eval();

        m_nees = (m_innovation.transpose() * S.colPivHouseholderQr().solve(m_innovation))(0, 0);
    }

    const state_vector_t &state() const
    {
        return m_x;
    }

    const cov_matrix_t &covariance() const
    {
        return m_P;
    }

    const output_vector_t &innovation() const
    {
        return m_innovation;
    }

    Scalar nees() const
    {
        return m_nees;
    }

private:
    Scalar m_eps;
    Scalar m_nees{0};
    Dynamics m_dynamics;
    cov_matrix_t m_P;
    cov_matrix_t m_Q;
    state_vector_t m_x;
    meas_cov_matrix_t m_R;
    Measurement m_measurement;
    output_vector_t m_innovation;
};

// CTAD deduction guide
template <typename Dynamics, typename Measurement, typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
ekf(Dynamics, Measurement, ekf_config<Scalar, NX, NU, NY>)
    -> ekf<Scalar, NX, NU, NY, Dynamics, Measurement>;

namespace detail {

struct ekf_sa_dynamics
{
    Vector<double, 2> operator()(const Vector<double, 2> &, const Vector<double, 1> &) const
    {
        return Vector<double, 2>::Zero();
    }
};

struct ekf_sa_measurement
{
    Vector<double, 1> operator()(const Vector<double, 2> &) const
    {
        return Vector<double, 1>::Zero();
    }
};

}

static_assert(ObserverPolicy<ekf<double, 2, 1, 1, detail::ekf_sa_dynamics, detail::ekf_sa_measurement>>);
static_assert(CovarianceObserver<ekf<double, 2, 1, 1, detail::ekf_sa_dynamics, detail::ekf_sa_measurement>>);

}

#endif
