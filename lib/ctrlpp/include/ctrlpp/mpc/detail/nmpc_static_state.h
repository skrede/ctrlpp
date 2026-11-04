#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_STATIC_STATE_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_STATIC_STATE_H

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/mpc/nlp_types.h"
#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/diagnostics.h"
#include "ctrlpp/mpc/nmpc_config.h"
#include "ctrlpp/mpc/nlp_formulation.h"

#include <Eigen/Dense>

#include <memory>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NH, typename Solver, typename Dynamics, std::size_t NC, std::size_t NTC>
class nmpc_static_state
{
public:
    static constexpr std::int32_t horizon = static_cast<std::int32_t>(NH);
    static constexpr int problem_dimension = static_cast<int>((NH + 1) * NX + NH * NU);

    using problem_type = nlp_problem_static<Scalar, problem_dimension>;

protected:
    Dynamics m_dynamics;
    nmpc_config<Scalar, NX, NU, NC, NTC> m_config;

    std::shared_ptr<nmpc_formulation_state<Scalar, NX, NU>> m_state;
    std::unique_ptr<problem_type> m_problem;
    Solver m_solver;

    Eigen::VectorX<Scalar> m_warm_z;
    Eigen::VectorX<Scalar> m_last_solution;
    Eigen::VectorX<Scalar> m_candidate;
    nlp_update<Scalar> m_update;
    mpc_diagnostics<Scalar> m_last_diagnostics;
    Vector<Scalar, NU> m_u_prev;
    bool m_setup_failed;
    bool m_has_solution;

    nmpc_static_state(Dynamics&& dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config)
        : m_dynamics{std::move(dynamics)}
        , m_config{config}
        , m_state{std::make_shared<nmpc_formulation_state<Scalar, NX, NU>>()}
        , m_solver{}
        , m_update{}
        , m_last_diagnostics{}
        , m_u_prev{Vector<Scalar, NU>::Zero()}
        , m_setup_failed{false}
        , m_has_solution{false}
    {
        initialize();
    }

    nmpc_static_state(Dynamics&& dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config, Solver&& solver)
        : m_dynamics{std::move(dynamics)}
        , m_config{config}
        , m_state{std::make_shared<nmpc_formulation_state<Scalar, NX, NU>>()}
        , m_solver{std::move(solver)}
        , m_update{}
        , m_last_diagnostics{}
        , m_u_prev{Vector<Scalar, NU>::Zero()}
        , m_setup_failed{false}
        , m_has_solution{false}
    {
        initialize();
    }

    // A defaulted move is correct: the problem sits behind a unique_ptr, so the reference the solver's bridge holds
    // to it survives the move.
    nmpc_static_state(nmpc_static_state&&) = default;
    nmpc_static_state& operator=(nmpc_static_state&&) = default;
    ~nmpc_static_state() = default;

    // A copy is an independent fork: it deep-copies the formulation state, poses its own problem against that copy
    // and sets its forked solver up against it.
    nmpc_static_state(const nmpc_static_state& other)
        : m_dynamics{other.m_dynamics}
        , m_config{other.m_config}
        , m_state{std::make_shared<nmpc_formulation_state<Scalar, NX, NU>>(*other.m_state)}
        , m_solver{other.m_solver}
        , m_warm_z{other.m_warm_z}
        , m_last_solution{other.m_last_solution}
        , m_candidate{other.m_candidate}
        , m_update{other.m_update}
        , m_last_diagnostics{other.m_last_diagnostics}
        , m_u_prev{other.m_u_prev}
        , m_setup_failed{false}
        , m_has_solution{other.m_has_solution}
    {
        build_problem_and_setup();
    }

    nmpc_static_state& operator=(const nmpc_static_state& other)
    {
        nmpc_static_state tmp{other};
        *this = std::move(tmp);
        return *this;
    }

    // Construction sizes the reused solve buffers once, so the steady-state assignments take the same-size path and
    // the solve allocates nothing.
    void initialize()
    {
        m_state->x_ref.resize(static_cast<std::size_t>(horizon) + 1, Vector<Scalar, NX>::Zero());
        m_warm_z = Eigen::VectorX<Scalar>::Zero(problem_dimension);
        m_update.x0 = Eigen::VectorX<Scalar>::Zero(problem_dimension);
        m_last_solution = Eigen::VectorX<Scalar>::Zero(problem_dimension);
        m_candidate = Eigen::VectorX<Scalar>::Zero(problem_dimension);
    }

    expected<void, nlp_formulation_error> pose_problem()
    {
        auto problem = detail::build_nmpc_problem_static<Scalar, NX, NU, NH, NC, NTC>(m_dynamics, m_config, m_state);
        if(!problem.has_value())
        {
            m_problem = std::make_unique<problem_type>();
            return unexpected(problem.error());
        }
        m_problem = std::make_unique<problem_type>(*std::move(problem));
        return {};
    }

    // A constructor cannot return an error, so a copy's refused re-pose or re-setup is the flag that turns its every
    // solve into solver_error::setup_incomplete.
    void build_problem_and_setup() { m_setup_failed = !pose_problem().has_value() || !detail::setup_nlp_solver(m_solver, *m_problem).has_value(); }
};

}

#endif
