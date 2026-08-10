#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_CT_CT_NLOC_ARM_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_CT_CT_NLOC_ARM_H

// The installed package advertises HPIPM in its interface definitions, but the
// interior-point backend that macro compiles calls d_ocp_qp_dim_set_all and
// d_ocp_qp_set_all with more arrays than the installed hpipm headers declare, so
// the umbrella does not compile under it. The backend selected below is the
// Gauss-Newton Riccati one, which the macro does not gate.
#undef HPIPM

#include "nloc_problem.h"

#include <cassert>  // must precede ct_optcon includes; DynamicRiccatiEquation.hpp uses assert() without <cassert>
#include <ct/optcon/optcon.h>

#include <Eigen/Dense>

#include <memory>
#include <cstdint>

namespace ctrlpp::bench
{

// A budget, not a tolerance: the achieved-cost row says whether it sufficed.
constexpr int32_t nloc_iteration_budget = 10;

// ct::core::LTISystem declares raw Eigen returns against a base whose returns
// are ct's own matrix wrappers, which is not a valid covariant override and
// fails to instantiate on every installed compiler, so the system is derived
// from the base with the base's own typedefs.
class linear_oscillator final : public ct::core::LinearSystem<oscillator_state_dim, oscillator_input_dim>
{
public:
    using base = ct::core::LinearSystem<oscillator_state_dim, oscillator_input_dim>;

    linear_oscillator()
        : base(), m_A{oscillator_dynamics_matrix()}, m_B{oscillator_input_matrix()}
    {
    }

    linear_oscillator* clone() const override
    {
        return new linear_oscillator(*this);
    }

    const typename base::state_matrix_t& getDerivativeState(const typename base::state_vector_t&,
                                                            const typename base::control_vector_t&,
                                                            const typename base::time_t) override
    {
        return m_A;
    }

    const typename base::state_control_matrix_t& getDerivativeControl(const typename base::state_vector_t&,
                                                                      const typename base::control_vector_t&,
                                                                      const typename base::time_t) override
    {
        return m_B;
    }

private:
    typename base::state_matrix_t         m_A;
    typename base::state_control_matrix_t m_B;
};

using nloc_cost = ct::optcon::CostFunctionQuadraticSimple<oscillator_state_dim, oscillator_input_dim>;
using nloc_problem_type = ct::optcon::ContinuousOptConProblem<oscillator_state_dim, oscillator_input_dim>;
using nloc_solver = ct::optcon::NLOptConSolver<oscillator_state_dim, oscillator_input_dim>;

inline std::shared_ptr<nloc_cost> make_nloc_cost()
{
    const ct::core::StateVector<oscillator_state_dim> x_nominal = ct::core::StateVector<oscillator_state_dim>::Zero();
    const ct::core::ControlVector<oscillator_input_dim> u_nominal =
        ct::core::ControlVector<oscillator_input_dim>::Zero();
    return std::make_shared<nloc_cost>(oscillator_state_weight::Identity(), oscillator_input_weight::Identity(),
                                       x_nominal, u_nominal, x_nominal, oscillator_state_weight::Identity());
}

inline ct::optcon::NLOptConSettings make_nloc_settings()
{
    ct::optcon::NLOptConSettings settings;
    settings.dt = oscillator_sample_period;
    settings.integrator = ct::core::IntegrationType::EULERCT;
    settings.discretization = ct::optcon::NLOptConSettings::APPROXIMATION::FORWARD_EULER;
    settings.max_iterations = nloc_iteration_budget;
    settings.nThreads = 1;
    settings.nlocp_algorithm = ct::optcon::NLOptConSettings::NLOCP_ALGORITHM::GNMS;
    settings.lqocp_solver = ct::optcon::NLOptConSettings::LQOCP_SOLVER::GNRICCATI_SOLVER;
    settings.printSummary = false;
    return settings;
}

inline nloc_solver::Policy_t make_nloc_guess()
{
    const std::size_t steps = static_cast<std::size_t>(oscillator_horizon);
    ct::core::StateVectorArray<oscillator_state_dim> states(steps + 1, oscillator_initial_state());
    ct::core::ControlVectorArray<oscillator_input_dim> feedforward(
        steps, ct::core::ControlVector<oscillator_input_dim>::Zero());
    ct::core::FeedbackArray<oscillator_state_dim, oscillator_input_dim> feedback(
        steps, ct::core::FeedbackMatrix<oscillator_state_dim, oscillator_input_dim>::Zero());
    return nloc_solver::Policy_t(states, feedforward, feedback, oscillator_sample_period);
}

// One closed-loop run per call, controller construction included, so a repeated
// call is a repeat of the same cold work rather than a warm continuation.
class ct_nloc_arm
{
public:
    ct_nloc_arm()
        : m_system{std::make_shared<linear_oscillator>()}
        , m_cost{make_nloc_cost()}
        , m_solver{nloc_problem_type{oscillator_horizon * oscillator_sample_period, oscillator_initial_state(),
                                     m_system, m_cost, m_system},
                   make_nloc_settings()}
    {
        m_solver.setInitialGuess(make_nloc_guess());
    }

    loop_outcome run()
    {
        oscillator_state x = oscillator_initial_state();
        loop_outcome out{0.0, 0.0};
        for(int32_t k = 0; k < oscillator_sim_steps; ++k)
        {
            m_solver.changeInitialState(x);
            m_solver.solve();
            const oscillator_input u = m_solver.getSolution().uff()[0];
            if(k == 0)
                out.first_input = u[0];
            out.cost += oscillator_stage_cost(x, u);
            x = oscillator_step(x, u);
        }
        return out;
    }

private:
    std::shared_ptr<linear_oscillator> m_system;
    std::shared_ptr<nloc_cost>         m_cost;
    nloc_solver                        m_solver;
};

}

#endif
