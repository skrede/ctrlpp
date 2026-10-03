#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_CLOSED_LOOP_BOUND_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_CLOSED_LOOP_BOUND_H

#include "derived_tolerance.h"

#include <Eigen/Dense>

#include <array>
#include <cstddef>

namespace ctrlpp {

/// @brief Bound on how far a computed run of x_{k+1} = (A - B*K)*x_k may land
/// from the exact run of a reference closed loop after `Steps` steps.
///
/// The reference is A_cl = A - B*K on the given matrices, run in exact
/// arithmetic from x0. A computed run uses its own A^, B^ and K^, so it follows
/// A^_cl = A_cl + E, and each step adds a rounding r_j. Telescoping gives
///   x^_k - x_k = sum_{j<k} A^_cl^(k-1-j) * (E*x_j + r_j),
///   A^_cl^m    = A_cl^m + sum_{i<m} A^_cl^(m-1-i) * E * A_cl^i,
/// so with c_m = |A_cl^m| the powers of the computed loop are bounded by
/// p_m = c_m + e * sum_{i<m} p_(m-1-i) * c_i, where e bounds |E|. The departure
/// is then bounded by d_k = sum_{j<k} p_(k-1-j) * (e*|x_j| + rho*eps*g*(|x_j| + d_j)),
/// where rho counts a step's roundings and g bounds the operands they land on.
/// No term is linearized away and the reference's decay enters through c_m and
/// |x_j|, which is what keeps the bound proportional to the decayed state.
///
/// Norms are infinity norms throughout. The recurrences are evaluated in
/// double; their own relative error is of order Steps^2 times double epsilon,
/// which this bound does not charge.
///
/// @cite higham2002 -- Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., 2002, Ch. 3 (componentwise rounding of products and sums)
template<int Nx, int Nu, std::size_t Steps>
class closed_loop_bound
{
public:
    using state_matrix = Eigen::Matrix<double, Nx, Nx>;
    using input_matrix = Eigen::Matrix<double, Nx, Nu>;
    using gain_matrix  = Eigen::Matrix<double, Nu, Nx>;
    using state_vector = Eigen::Vector<double, Nx>;

    closed_loop_bound(const state_matrix &A, const input_matrix &B, const gain_matrix &K, const state_vector &x0)
            : a_norm_(A.cwiseAbs().rowwise().sum().maxCoeff())
            , b_norm_(B.cwiseAbs().rowwise().sum().maxCoeff())
            , k_norm_(K.cwiseAbs().rowwise().sum().maxCoeff())
            , departure_{}
            , power_norm_{}
            , state_norm_{}
            , perturbed_power_{}
    {
        const state_matrix closed = A - B * K;
        state_matrix power        = state_matrix::Identity();
        state_vector state        = x0;
        for(std::size_t m = 0; m <= Steps; ++m)
        {
            power_norm_[m] = power.cwiseAbs().rowwise().sum().maxCoeff();
            state_norm_[m] = state.cwiseAbs().maxCoeff();
            power          = closed * power;
            state          = closed * state;
        }
    }

    /// `eps` is the epsilon of the run's arithmetic, `plant_roundings` bounds
    /// the roundings that formed each entry of its A^ and B^ from the
    /// reference's, and `gain_departure` bounds |K^ - K|.
    double final_state_departure(double eps, double plant_roundings, double gain_departure)
    {
        const double gain   = k_norm_ + gain_departure;
        const double plant  = plant_roundings * eps;
        const double e      = plant * a_norm_ + plant * b_norm_ * gain + b_norm_ * gain_departure;
        const double g      = (1.0 + plant) * (a_norm_ + b_norm_ * gain);
        const double charge = loop_roundings_per_step(Nx, Nu) * eps * g;
        fill_perturbed_power(e);
        departure_[0] = 0.0;
        for(std::size_t k = 1; k <= Steps; ++k)
        {
            double sum = 0.0;
            for(std::size_t j = 0; j < k; ++j)
                sum += perturbed_power_[k - 1 - j] * (e * state_norm_[j] + charge * (state_norm_[j] + departure_[j]));
            departure_[k] = sum;
        }
        return departure_[Steps];
    }

private:
    double a_norm_;
    double b_norm_;
    double k_norm_;
    std::array<double, Steps + 1> departure_;
    std::array<double, Steps + 1> power_norm_;
    std::array<double, Steps + 1> state_norm_;
    std::array<double, Steps + 1> perturbed_power_;

    void fill_perturbed_power(double e)
    {
        perturbed_power_[0] = 1.0;
        for(std::size_t m = 1; m <= Steps; ++m)
        {
            double sum = 0.0;
            for(std::size_t i = 0; i < m; ++i)
                sum += perturbed_power_[m - 1 - i] * power_norm_[i];
            perturbed_power_[m] = power_norm_[m] + e * sum;
        }
    }
};

}

#endif
