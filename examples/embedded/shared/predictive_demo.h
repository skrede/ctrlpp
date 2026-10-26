#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_PREDICTIVE_DEMO_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_PREDICTIVE_DEMO_H

#include "golden_reference.h"

#include "ctrlpp/nmpc.h"
#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/nmpc_config.h"
#include "ctrlpp/mpc/argmin_solver.h"

#include <type_traits>

namespace ctrlpp {

inline const char *describe(solver_error e)
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

template<class Scalar>
struct predictive_plant
{
    Vector<Scalar, kPredictiveNx> operator()(const Vector<Scalar, kPredictiveNx> &x, const Vector<Scalar, kPredictiveNu> &u) const
    {
        const Scalar dt = static_cast<Scalar>(kPredictiveDt);
        return Vector<Scalar, kPredictiveNx>{x[0] + dt * x[1], x[1] + dt * u[0]};
    }
};

// The configuration nmpc_static_nomalloc_test pins and the embedded compile
// witness builds: the constraint bound is the equality count, which keeps the
// solver's per-call multiplier storage inline. The backend's NW-SQP policy fixes
// its scalar to double, so the family has no single-precision form.
template<class Scalar>
class predictive_demo
{
    static_assert(std::is_same_v<Scalar, double>, "the backend's NW-SQP policy computes in double only");

public:
    using plant_type      = predictive_plant<Scalar>;
    using state_type      = Vector<Scalar, kPredictiveNx>;
    using input_type      = Vector<Scalar, kPredictiveNu>;
    using solver_type     = argmin_solver<Scalar, argmin_nw_sqp, true, kPredictiveNv, kPredictiveMaxM>;
    using controller_type = nmpc_static<Scalar, kPredictiveNx, kPredictiveNu, kPredictiveNh, solver_type, plant_type>;

    // The controller carries its solver workspace inline, about 42 KB here. The
    // factory's result initializes the member directly, so the controller is
    // built once, in place, and a refused configuration is a value the demo
    // reports rather than a controller it runs.
    predictive_demo()
            : controller_(controller_type::create(plant_type{}, config()))
            , x_(Scalar{1}, Scalar{0})
            , u_(input_type::Zero())
    {
    }

    static nmpc_config<Scalar, kPredictiveNx, kPredictiveNu> config()
    {
        nmpc_config<Scalar, kPredictiveNx, kPredictiveNu> settings;
        settings.horizon = static_cast<int>(kPredictiveNh);
        settings.Q       = static_cast<Scalar>(kPredictiveStateWeight) * Matrix<Scalar, kPredictiveNx, kPredictiveNx>::Identity();
        settings.R       = static_cast<Scalar>(kPredictiveInputWeight) * Matrix<Scalar, kPredictiveNu, kPredictiveNu>::Identity();
        return settings;
    }

    // One solve followed by the plant step: the block the host allocation test
    // arms, so the board's figure and the host's cover the same body.
    ctrlpp::expected<void, solver_error> advance()
    {
        if(!constructed())
            return ctrlpp::unexpected(solver_error::setup_incomplete);
        const auto solved = controller_->solve(x_);
        if(!solved.has_value())
            return ctrlpp::unexpected(solved.error());
        u_ = solved->input;
        x_ = plant_type{}(x_, u_);
        return {};
    }

    bool constructed() const
    {
        return controller_.has_value();
    }

    // Valid only when constructed().
    const controller_type &controller() const
    {
        return *controller_;
    }

    const state_type &state() const
    {
        return x_;
    }

    const input_type &input() const
    {
        return u_;
    }

private:
    ctrlpp::expected<controller_type, typename controller_type::construction_error> controller_;
    state_type x_;
    input_type u_;
};

}

#endif
