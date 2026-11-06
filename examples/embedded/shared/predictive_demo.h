#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_PREDICTIVE_DEMO_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_PREDICTIVE_DEMO_H

#include "golden_reference.h"

#include "ctrlpp/nmpc.h"
#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/nmpc_config.h"
#include "ctrlpp/mpc/argmin_solver.h"

#include <cstdint>
#include <type_traits>

namespace ctrlpp
{

inline const char* describe(solver_error e)
{
    switch(e)
    {
        case solver_error::infeasible:
            return "the problem is infeasible";
        case solver_error::invalid_problem:
            return "the solver rejected the problem";
        case solver_error::setup_incomplete:
            return "the formulation was refused at setup";
        case solver_error::invalid_backend_result:
            return "the backend returned a malformed primal";
    }
    return "unknown";
}

template <class Scalar>
struct predictive_plant
{
    Vector<Scalar, golden::predictive_nx> operator()(const Vector<Scalar, golden::predictive_nx>& x, const Vector<Scalar, golden::predictive_nu>& u) const
    {
        const Scalar dt = static_cast<Scalar>(golden::predictive_dt);
        return Vector<Scalar, golden::predictive_nx>{x[0] + dt * x[1], x[1] + dt * u[0]};
    }
};

// The configuration the host no-allocation test pins and the embedded compile witness builds: bounding
// the constraints at the equality count keeps the solver's per-call multiplier storage inline.
template <class Scalar>
class predictive_demo
{
    static_assert(std::is_same_v<Scalar, double>, "the backend's NW-SQP policy computes in double only");

public:
    using plant_type = predictive_plant<Scalar>;
    using state_type = Vector<Scalar, golden::predictive_nx>;
    using input_type = Vector<Scalar, golden::predictive_nu>;
    using solver_type = argmin_solver<Scalar, argmin_nw_sqp, true, golden::predictive_nv, golden::predictive_max_m>;
    using controller_type = nmpc_static<Scalar, golden::predictive_nx, golden::predictive_nu, golden::predictive_nh, solver_type, plant_type>;

    // The controller carries its solver workspace inline, about 42 KB here, so the factory's result
    // initializes the member in place; a refused configuration is a value the demo reports, not one it runs.
    predictive_demo()
        : m_controller(controller_type::create(plant_type{}, config()))
        , m_x(Scalar{1}, Scalar{0})
        , m_u(input_type::Zero())
    {
    }

    static nmpc_config<Scalar, golden::predictive_nx, golden::predictive_nu> config()
    {
        nmpc_config<Scalar, golden::predictive_nx, golden::predictive_nu> settings;
        settings.horizon = static_cast<std::int32_t>(golden::predictive_nh);
        settings.Q = static_cast<Scalar>(golden::predictive_state_weight) * Matrix<Scalar, golden::predictive_nx, golden::predictive_nx>::Identity();
        settings.R = static_cast<Scalar>(golden::predictive_input_weight) * Matrix<Scalar, golden::predictive_nu, golden::predictive_nu>::Identity();
        return settings;
    }

    // One solve and the plant step, the block the host allocation test arms, so the board's figure and the
    // host's cover the same body.
    ctrlpp::expected<void, solver_error> advance()
    {
        if(!constructed())
            return ctrlpp::unexpected(solver_error::setup_incomplete);
        const auto solved = m_controller->solve(m_x);
        if(!solved.has_value())
            return ctrlpp::unexpected(solved.error());
        m_u = solved->input;
        m_x = plant_type{}(m_x, m_u);
        return {};
    }

    bool constructed() const { return m_controller.has_value(); }

    // Valid only when constructed().
    const controller_type& controller() const { return *m_controller; }

    const state_type& state() const { return m_x; }

    const input_type& input() const { return m_u; }

private:
    ctrlpp::expected<controller_type, typename controller_type::construction_error> m_controller;
    state_type m_x;
    input_type m_u;
};

}

#endif
