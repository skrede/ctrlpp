// Embedded-clean compile witness for the predictive surface: instantiates
// nmpc_static over the fixed-dimension argmin solver at the configuration the
// allocation test pins, and exercises one construction and one solve. Every
// fallible result is handled through has_value() plus operator*, so the unit
// stays free of the exception machinery and compiles under -fno-exceptions with
// CTRLPP_NO_EXCEPTIONS defined. It needs the optional nonlinear-programming
// backend, so it is registered only when that backend is on.
//
// The scalar is double because the backend's NW-SQP policy fixes its scalar type
// to double; the same instantiation at float does not compile.

#include "ctrlpp/nmpc.h"
#include "ctrlpp/types.h"

#include "ctrlpp/mpc/argmin_solver.h"

#include <cstddef>

namespace
{

using scalar = double;

constexpr std::size_t nx = 2;
constexpr std::size_t nu = 1;
constexpr std::size_t nh = 5;

constexpr int decision_dimension = static_cast<int>((nh + 1) * nx + nh * nu);

// Bound to the equality count -- dynamics continuity plus the initial state --
// so the solver's per-call multiplier storage stays inline.
constexpr int constraint_dimension = static_cast<int>(nx * (nh + 1));

constexpr scalar dt = 0.1;

struct double_integrator
{
    auto operator()(const ctrlpp::Vector<scalar, nx>& x, const ctrlpp::Vector<scalar, nu>& u) const
        -> ctrlpp::Vector<scalar, nx>
    {
        ctrlpp::Vector<scalar, nx> next;
        next[0] = x[0] + dt * x[1];
        next[1] = x[1] + dt * u[0];
        return next;
    }
};

using solver_type = ctrlpp::argmin_solver<scalar, ctrlpp::argmin_nw_sqp, true, decision_dimension, constraint_dimension>;
using controller_type = ctrlpp::nmpc_static<scalar, nx, nu, nh, solver_type, double_integrator>;

static_assert(controller_type::problem_dimension == decision_dimension);

}

int main()
{
    ctrlpp::nmpc_config<scalar, nx, nu> config;
    config.horizon = static_cast<int>(nh);

    controller_type controller{double_integrator{}, config};

    ctrlpp::Vector<scalar, nx> x0 = ctrlpp::Vector<scalar, nx>::Zero();
    x0[0] = scalar{1};

    const auto u = controller.solve(x0);
    return (u.has_value() && (*u).input.allFinite()) ? 0 : 1;
}
