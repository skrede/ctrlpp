#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_CT_CT_LQR_ARM_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_CT_CT_LQR_ARM_H

#include "riccati_problem.h"

#include <cassert>  // must precede ct_optcon includes; DynamicRiccatiEquation.hpp uses assert() without <cassert>
#include <ct/core/types/StateVector.h>
#include <ct/core/types/ControlVector.h>
#include <ct/optcon/lqr/riccati/CARE.hpp>
#include <ct/optcon/lqr/riccati/CARE-impl.hpp>
#include <ct/optcon/lqr/LQR.hpp>
#include <ct/optcon/lqr/LQR-impl.hpp>

#include <Eigen/Dense>

#include <cstddef>

namespace ctrlpp::bench
{

// ct's LQR writes its gain into a caller-owned matrix, so the arm holds that
// matrix rather than returning one: a return by value inside the timed region
// would charge the ct arm for a copy the ctrlpp arm never makes.
template <std::size_t NX, std::size_t NU>
class ct_lqr_arm
{
public:
    explicit ct_lqr_arm(const riccati_plant<NX, NU>& plant)
        : m_lqr{}, m_K{}, m_B{plant.B}, m_A{plant.A}, m_Q{plant.Q}, m_R{plant.R}
    {
    }

    void solve()
    {
        m_lqr.compute(m_Q, m_R, m_A, m_B, m_K);
    }

    const Eigen::Matrix<double, int(NU), int(NX)>& gain() const
    {
        return m_K;
    }

private:
    ct::optcon::LQR<NX, NU>                            m_lqr;
    Eigen::Matrix<double, int(NU), int(NX)>            m_K;
    Eigen::Matrix<double, int(NX), int(NU)>            m_B;
    typename ct::optcon::LQR<NX, NU>::state_matrix_t   m_A;
    typename ct::optcon::LQR<NX, NU>::state_matrix_t   m_Q;
    typename ct::optcon::LQR<NX, NU>::control_matrix_t m_R;
};

}

#endif
