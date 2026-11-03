#ifndef HPP_GUARD_CTRLPP_ESTIMATION_LUENBERGER_H
#define HPP_GUARD_CTRLPP_ESTIMATION_LUENBERGER_H

// @cite luenberger1971 -- Luenberger, "An Introduction to Observers", 1971

#include "ctrlpp/estimation/observer_policy.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/model/state_space.h"

#include <cstddef>
#include <utility>

namespace ctrlpp
{

enum class luenberger_update_error
{
    non_finite_state,
    non_finite_measurement,
    non_finite_result,
};

enum class luenberger_health
{
    ok,
    non_finite_estimate,
};

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
class luenberger_observer
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NY > 0, "Output dimension NY must be positive");

    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);
    static constexpr int ny = static_cast<int>(NY);

public:
    using observer_tag = struct luenberger_tag;
    using state_vector_t = Eigen::Matrix<Scalar, nx, 1>;
    using input_vector_t = Eigen::Matrix<Scalar, nu, 1>;
    using output_vector_t = Eigen::Matrix<Scalar, ny, 1>;
    using gain_matrix_t = Eigen::Matrix<Scalar, nx, ny>;
    using system_t = discrete_state_space<Scalar, NX, NU, NY>;

    luenberger_observer(system_t sys, gain_matrix_t L, state_vector_t x0)
        : m_sys{std::move(sys)}
        , m_L{std::move(L)}
        , m_x{std::move(x0)}
        , m_health{luenberger_health::ok}
    {
    }

    // Not fallible: its input is a command the plant already took. A state it poisons is
    // latched into health() by the next update.
    void predict(const input_vector_t& u) { m_x = (m_sys.A * m_x + m_sys.B * u).eval(); }

    ctrlpp::expected<void, luenberger_update_error> update(const output_vector_t& z)
    {
        if(const auto step = check_step(z); !step)
            return ctrlpp::unexpected(latch_health(step.error()));

        auto next_x = (m_x + m_L * (z - m_sys.C * m_x)).eval();
        if(!next_x.allFinite())
            return ctrlpp::unexpected(luenberger_update_error::non_finite_result);
        m_x = std::move(next_x);
        return {};
    }

    const state_vector_t& state() const { return m_x; }

    luenberger_health health() const { return m_health; }

    void set_gain(const gain_matrix_t& L) { m_L = L; }
    void set_model(system_t sys) { m_sys = std::move(sys); }

    void reset(const state_vector_t& x0)
    {
        m_x = x0;
        m_health = luenberger_health::ok;
    }

private:
    system_t m_sys;
    gain_matrix_t m_L;
    state_vector_t m_x;
    luenberger_health m_health;

    // Classifies without touching a member: the carried estimate first, the measurement last.
    ctrlpp::expected<void, luenberger_update_error> check_step(const output_vector_t& z) const
    {
        if(!m_x.allFinite())
            return ctrlpp::unexpected(luenberger_update_error::non_finite_state);
        if(!z.allFinite())
            return ctrlpp::unexpected(luenberger_update_error::non_finite_measurement);
        return {};
    }

    luenberger_update_error latch_health(luenberger_update_error fault)
    {
        if(fault != luenberger_update_error::non_finite_measurement)
            m_health = luenberger_health::non_finite_estimate;
        return fault;
    }
};

static_assert(ObserverPolicy<luenberger_observer<double, 2, 1, 1>>);
static_assert(!CovarianceObserver<luenberger_observer<double, 2, 1, 1>>);

}

#endif
