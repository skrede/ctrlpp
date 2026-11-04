#ifndef HPP_GUARD_CTRLPP_MHE_DETAIL_NMHE_TERMS_H
#define HPP_GUARD_CTRLPP_MHE_DETAIL_NMHE_TERMS_H

#include "ctrlpp/mhe/detail/nmhe_layout.h"

#include "ctrlpp/mhe/mhe_config.h"

#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <span>
#include <memory>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <functional>

namespace ctrlpp::detail
{

// Members copied from a const source are const, as a lambda's by-copy captures of them would be, so handing the
// closure to std::function copies them rather than moving them.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N, std::size_t NC, typename Dynamics, typename Measurement, typename State>
class nmhe_cost
{
    using layout_type = nmhe_layout<NX, N>;

    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nc = static_cast<int>(NC);
    static constexpr int Ni = static_cast<int>(N);

public:
    nmhe_cost(const Dynamics& dynamics, const Measurement& measurement, const nmhe_config<Scalar, NX, NU, NY, N, NC>& config, std::shared_ptr<State> state,
            const Matrix<Scalar, NX, NX>& Q_inv, const Matrix<Scalar, NY, NY>& R_inv, const layout_type& layout)
        : m_has_path_slack{layout.has_path_slack}
        , m_arrival_weight{config.arrival_cost_weight}
        , m_state{std::move(state)}
        , m_dynamics{dynamics}
        , m_measurement{measurement}
        , m_Q_inv{Q_inv}
        , m_R_inv{R_inv}
        , m_config{config}
    {
    }

    Scalar operator()(std::span<const Scalar> z) const
    {
        Scalar total{0};
        Eigen::Map<const Vector<Scalar, NX>> x0(z.data() + layout_type::x_offset);
        Vector<Scalar, NX> dx0 = x0 - m_state->arrival_state;
        total += Scalar{0.5} * dx0.dot((m_arrival_weight * m_state->arrival_P_inv) * dx0);

        for(std::int32_t k = 0; k < Ni; ++k)
        {
            Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + layout_type::x_offset + k * nx);
            Eigen::Map<const Vector<Scalar, NX>> xk1(z.data() + layout_type::x_offset + (k + 1) * nx);
            Vector<Scalar, NX> wk = xk1 - m_dynamics(xk, m_state->u_window[static_cast<std::size_t>(k)]);
            total += Scalar{0.5} * wk.dot(m_Q_inv * wk);
        }
        for(std::int32_t k = 0; k <= Ni; ++k)
        {
            Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + layout_type::x_offset + k * nx);
            Vector<Scalar, NY> vk = m_state->z_window[static_cast<std::size_t>(k)] - m_measurement(xk);
            total += Scalar{0.5} * vk.dot(m_R_inv * vk);
        }
        if(m_has_path_slack)
            for(std::int32_t k = 0; k <= Ni; ++k)
                total += m_config.path_penalty.dot(Eigen::Map<const Vector<Scalar, NC>>(z.data() + layout_type::path_slack_offset + k * nc));

        return total;
    }

private:
    const bool m_has_path_slack;
    const Scalar m_arrival_weight;
    std::shared_ptr<State> m_state;
    const Dynamics m_dynamics;
    const Measurement m_measurement;
    const Matrix<Scalar, NX, NX> m_Q_inv;
    const Matrix<Scalar, NY, NY> m_R_inv;
    const nmhe_config<Scalar, NX, NU, NY, N, NC> m_config;
};

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N, std::size_t NC, typename Measurement, typename State>
class nmhe_constraints
{
    using layout_type = nmhe_layout<NX, N>;
    using path_constraint_type = std::function<Vector<Scalar, NC>(const Vector<Scalar, NX>&)>;

    static constexpr int nx = static_cast<int>(NX);
    static constexpr int ny = static_cast<int>(NY);
    static constexpr int nc = static_cast<int>(NC);
    static constexpr int Ni = static_cast<int>(N);

public:
    nmhe_constraints(const Measurement& measurement, const nmhe_config<Scalar, NX, NU, NY, N, NC>& config, std::shared_ptr<State> state, const layout_type& layout)
        : m_has_path_slack{layout.has_path_slack}
        , m_has_residual{layout.has_residual}
        , m_residual_con_start{layout.residual_con_start}
        , m_state{std::move(state)}
        , m_measurement{measurement}
        , m_config{config}
    {
    }

    void operator()(std::span<const Scalar> z, std::span<Scalar> c) const
    {
        if constexpr(NC > 0)
        {
            if(m_config.path_constraint_cb)
            {
                const auto& g = *m_config.path_constraint_cb;
                for(std::int32_t k = 0; k <= Ni; ++k)
                    set_path_stage(g, z, c, k);
            }
        }

        if(m_has_residual)
        {
            const auto& bound = *m_config.residual_bound;
            for(std::int32_t k = 0; k <= Ni; ++k)
                set_residual_stage(bound, z, c, k);
        }
    }

private:
    const bool m_has_path_slack;
    const bool m_has_residual;
    const Eigen::Index m_residual_con_start;
    std::shared_ptr<State> m_state;
    const Measurement m_measurement;
    const nmhe_config<Scalar, NX, NU, NY, N, NC> m_config;

    void set_path_stage(const path_constraint_type& g, std::span<const Scalar> z, std::span<Scalar> c, std::int32_t k) const
    {
        Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + layout_type::x_offset + k * nx);

        Vector<Scalar, NC> gk = g(xk);

        if(m_has_path_slack)
        {
            Eigen::Map<const Vector<Scalar, NC>> sk(z.data() + layout_type::path_slack_offset + k * nc);
            gk -= sk;
        }

        for(std::int32_t i = 0; i < nc; ++i)
            c[static_cast<std::size_t>(layout_type::path_con_start + k * nc + i)] = gk[i];
    }

    void set_residual_stage(const Vector<Scalar, NY>& bound, std::span<const Scalar> z, std::span<Scalar> c, std::int32_t k) const
    {
        Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + layout_type::x_offset + k * nx);

        Vector<Scalar, NY> residual = m_state->z_window[static_cast<std::size_t>(k)] - m_measurement(xk);

        // Without a path block the residual rows start at row zero, a constant the compiler then folds.
        const Eigen::Index row_start = NC > 0 ? m_residual_con_start : 0;
        for(std::int32_t i = 0; i < ny; ++i)
        {
            std::size_t base = static_cast<std::size_t>(row_start + k * ny * 2 + i * 2);
            c[base] = residual[i] - bound[i];
            c[base + 1] = -residual[i] - bound[i];
        }
    }
};

}

#endif
