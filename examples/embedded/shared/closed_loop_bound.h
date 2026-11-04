#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_CLOSED_LOOP_BOUND_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_CLOSED_LOOP_BOUND_H

#include "derived_tolerance.h"

#include <Eigen/Dense>

#include <array>
#include <cstddef>

namespace ctrlpp
{

// Bounds how far a computed run of x_{k+1} = (A - B*K)*x_k lands from the exact run of the reference closed
// loop A_cl = A - B*K, started at x0, after Steps steps. A computed run follows A^_cl = A_cl + E and adds a
// rounding r_j each step; telescoping gives x^_k - x_k = sum_{j<k} A^_cl^(k-1-j) * (E*x_j + r_j) and
// A^_cl^m = A_cl^m + sum_{i<m} A^_cl^(m-1-i) * E * A_cl^i, so with c_m = |A_cl^m| and e >= |E| the computed
// powers are bounded by p_m = c_m + e * sum_{i<m} p_(m-1-i) * c_i, and the departure by
// d_k = sum_{j<k} p_(k-1-j) * (e*|x_j| + rho*eps*g*(|x_j| + d_j)), where rho counts a step's roundings and g
// bounds the operands they land on. No term is linearized away, and the reference's decay enters through
// c_m and |x_j|, which keeps the bound proportional to the decayed state. Norms are infinity norms; the
// recurrences run in double, and their own relative error, of order Steps^2 double epsilons, is not charged.
// @cite higham2002 -- Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., 2002, Ch. 3 (componentwise rounding of products and sums)
template <int Nx, int Nu, std::size_t Steps>
class closed_loop_bound
{
public:
    using state_matrix = Eigen::Matrix<double, Nx, Nx>;
    using input_matrix = Eigen::Matrix<double, Nx, Nu>;
    using gain_matrix = Eigen::Matrix<double, Nu, Nx>;
    using state_vector = Eigen::Vector<double, Nx>;

    closed_loop_bound(const state_matrix& A, const input_matrix& B, const gain_matrix& K, const state_vector& x0)
        : m_a_norm(A.cwiseAbs().rowwise().sum().maxCoeff())
        , m_b_norm(B.cwiseAbs().rowwise().sum().maxCoeff())
        , m_k_norm(K.cwiseAbs().rowwise().sum().maxCoeff())
        , m_departure{}
        , m_power_norm{}
        , m_state_norm{}
        , m_perturbed_power{}
    {
        const state_matrix closed = A - B * K;
        state_matrix power = state_matrix::Identity();
        state_vector state = x0;
        for(std::size_t m = 0; m <= Steps; ++m)
        {
            m_power_norm[m] = power.cwiseAbs().rowwise().sum().maxCoeff();
            m_state_norm[m] = state.cwiseAbs().maxCoeff();
            power = closed * power;
            state = closed * state;
        }
    }

    // plant_roundings bounds the roundings that formed each entry of the run's A^ and B^ from the
    // reference's, and gain_departure bounds |K^ - K|.
    double final_state_departure(double eps, double plant_roundings, double gain_departure)
    {
        const double gain = m_k_norm + gain_departure;
        const double plant = plant_roundings * eps;
        const double e = plant * m_a_norm + plant * m_b_norm * gain + m_b_norm * gain_departure;
        const double g = (1.0 + plant) * (m_a_norm + m_b_norm * gain);
        const double charge = loop_roundings_per_step(Nx, Nu) * eps * g;
        fill_perturbed_power(e);
        m_departure[0] = 0.0;
        for(std::size_t k = 1; k <= Steps; ++k)
        {
            double sum = 0.0;
            for(std::size_t j = 0; j < k; ++j)
                sum += m_perturbed_power[k - 1 - j] * (e * m_state_norm[j] + charge * (m_state_norm[j] + m_departure[j]));
            m_departure[k] = sum;
        }
        return m_departure[Steps];
    }

    // The departure of every step of a run whose input departs from the reference law by at most
    // offset[j] + slope * |x_j| at step j, on the reference's own A, B and K.
    const std::array<double, Steps + 1>& input_departures(double eps, const std::array<double, Steps>& offset, double slope)
    {
        const double g = m_a_norm + m_b_norm * m_k_norm;
        const double charge = loop_roundings_per_step(Nx, Nu) * eps * g;
        fill_perturbed_power(0.0);
        m_departure[0] = 0.0;
        for(std::size_t k = 1; k <= Steps; ++k)
        {
            double sum = 0.0;
            for(std::size_t j = 0; j < k; ++j)
            {
                const double state = m_state_norm[j] + m_departure[j];
                sum += m_perturbed_power[k - 1 - j] * (charge * state + m_b_norm * (offset[j] + slope * state));
            }
            m_departure[k] = sum;
        }
        return m_departure;
    }

    const std::array<double, Steps + 1>& reference_state_norms() const { return m_state_norm; }

    double gain_norm() const { return m_k_norm; }

private:
    double m_a_norm;
    double m_b_norm;
    double m_k_norm;
    std::array<double, Steps + 1> m_departure;
    std::array<double, Steps + 1> m_power_norm;
    std::array<double, Steps + 1> m_state_norm;
    std::array<double, Steps + 1> m_perturbed_power;

    void fill_perturbed_power(double e)
    {
        m_perturbed_power[0] = 1.0;
        for(std::size_t m = 1; m <= Steps; ++m)
        {
            double sum = 0.0;
            for(std::size_t i = 0; i < m; ++i)
                sum += m_perturbed_power[m - 1 - i] * m_power_norm[i];
            m_perturbed_power[m] = m_power_norm[m] + e * sum;
        }
    }
};

}

#endif
