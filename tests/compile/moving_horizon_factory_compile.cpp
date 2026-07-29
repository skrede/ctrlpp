#include "ctrlpp/mhe.h"
#include "ctrlpp/nmhe.h"

#include <cstddef>
#include <utility>
#include <type_traits>

namespace
{

using moving_horizon_error = ctrlpp::moving_horizon_configuration_error;
using matrix_type = ctrlpp::Matrix<double, 2, 2>;
using config_type = ctrlpp::mhe_config<double, 2, 1, 1, 4>;

using inverse_result = decltype(ctrlpp::detail::finite_full_piv_inverse(
    std::declval<const matrix_type&>(),
    moving_horizon_error::non_invertible_process_noise));
using validation_result =
    decltype(ctrlpp::detail::validate_moving_horizon_options(
        std::declval<const config_type&>()));

static_assert(std::is_same_v<
              inverse_result,
              ctrlpp::expected<matrix_type, moving_horizon_error>>);
static_assert(std::is_same_v<
              validation_result,
              ctrlpp::expected<void, moving_horizon_error>>);

enum class setup_error
{
    failed,
};

struct qp_solver_stub
{
    using scalar_type = double;

    auto setup(const ctrlpp::qp_problem<double>&)
        -> ctrlpp::expected<void, setup_error>
    {
        return {};
    }

    auto solve(const ctrlpp::qp_update<double>&)
        -> ctrlpp::qp_result<double>
    {
        return {};
    }
};

struct nlp_solver_stub
{
    using scalar_type = double;

    auto setup(const ctrlpp::nlp_problem<double>&)
        -> ctrlpp::expected<void, setup_error>
    {
        return {};
    }

    auto solve(const ctrlpp::nlp_update<double>&)
        -> ctrlpp::nlp_result<double>
    {
        return {};
    }
};

struct constant_dynamics
{
    auto operator()(const ctrlpp::Vector<double, 2>& x,
                    const ctrlpp::Vector<double, 1>&) const
        -> ctrlpp::Vector<double, 2>
    {
        return x;
    }
};

struct position_measurement
{
    auto operator()(const ctrlpp::Vector<double, 2>& x) const
        -> ctrlpp::Vector<double, 1>
    {
        return ctrlpp::Vector<double, 1>{x[0]};
    }
};

}

int main()
{
    constexpr std::size_t nx = 2;
    constexpr std::size_t nu = 1;
    constexpr std::size_t ny = 1;
    constexpr std::size_t horizon = 4;

    using mhe_type =
        ctrlpp::mhe<double, nx, nu, ny, horizon, qp_solver_stub,
                    constant_dynamics, position_measurement>;
    auto mhe_result = mhe_type::create(
        constant_dynamics{}, position_measurement{},
        ctrlpp::mhe_config<double, nx, nu, ny, horizon>{});
    if(!mhe_result)
        return 1;
    auto linear_estimator = std::move(*mhe_result);

    using nmhe_type =
        ctrlpp::nmhe<double, nx, nu, ny, horizon, nlp_solver_stub,
                     constant_dynamics, position_measurement>;
    auto nmhe_result = nmhe_type::create(
        constant_dynamics{}, position_measurement{},
        ctrlpp::nmhe_config<double, nx, nu, ny, horizon>{});
    if(!nmhe_result)
        return 1;
    auto nonlinear_estimator = std::move(*nmhe_result);

    linear_estimator.predict(ctrlpp::Vector<double, nu>::Zero());
    nonlinear_estimator.predict(ctrlpp::Vector<double, nu>::Zero());
    return 0;
}
