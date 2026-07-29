#ifndef HPP_GUARD_CTRLPP_ESTIMATION_LUENBERGER_H
#define HPP_GUARD_CTRLPP_ESTIMATION_LUENBERGER_H

/// @brief Discrete-time Luenberger state observer with fixed gain.
///
/// @cite luenberger1971 -- Luenberger, "An Introduction to Observers", 1971

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/model/state_space.h"

#include "ctrlpp/estimation/observer_policy.h"

#include <cstddef>
#include <utility>

namespace ctrlpp
{

/// @brief Structured failure modes of a `luenberger_observer` measurement
/// update.
///
/// Each enumerator is an exact domain condition, not a tuning preference: a
/// non-finite operand makes every downstream product non-finite, so the step
/// cannot produce an estimate at all. The observer carries no covariance, so
/// there are two causes rather than the three the covariance filters have.
///
///  * non_finite_state       : the carried state estimate is already non-finite
///                             when the step begins. It is reported ahead of the
///                             measurement because the fault is upstream of it:
///                             a caller told "your measurement is bad" would
///                             replace a working sensor while the real fault
///                             sits in the prediction that poisoned the state.
///  * non_finite_measurement : the supplied measurement vector has a non-finite
///                             component. The correction is the single fused
///                             expression x + L*(z - C*x), so every state
///                             component whose gain row is nonzero becomes
///                             non-finite, and the observer carries that state
///                             forward with no path back.
enum class luenberger_update_error
{
    non_finite_state,
    non_finite_measurement,
    non_finite_result,
};

/// @brief Persistent state-health status of a `luenberger_observer`.
///
/// A per-call result cannot answer whether the carried estimate is still
/// degraded from a step several samples ago, because that question outlives the
/// call. The status latches until `reset` seeds a fresh estimate.
///
///  * ok                  : every step so far began from a finite estimate.
///  * non_finite_estimate : a step found the carried state already non-finite.
///                          `predict` does not reject its input, so this is how
///                          a poisoned control vector or a non-finite model
///                          becomes visible. A rejected measurement does NOT set
///                          it: the rejection mutates nothing, so it leaves the
///                          observer healthy.
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

    luenberger_observer(system_t sys, gain_matrix_t L, state_vector_t x0) : m_sys{std::move(sys)}, m_L{std::move(L)}, m_x{std::move(x0)} {}

    void predict(const input_vector_t& u) { m_x = (m_sys.A * m_x + m_sys.B * u).eval(); }

    /// @brief Correct the carried state with a measurement.
    ///
    /// The step is rejected before the state is assigned when the carried
    /// estimate or the measurement is non-finite, so a rejected step leaves the
    /// state bitwise unchanged and the caller may retry with the next sample.
    ///
    /// `predict` is deliberately not fallible: its input is a control vector the
    /// caller already commanded and owns, and rejecting it would leave the
    /// observer with no propagation for a step the plant did take. A prediction
    /// that poisons the state is instead reported by `health()`, which the next
    /// update latches.
    auto update(const output_vector_t& z) -> ctrlpp::expected<void, luenberger_update_error>
    {
        if(const auto step = check_step(z); !step)
            return ctrlpp::unexpected(latch_health(step.error()));

        auto next_x = (m_x + m_L * (z - m_sys.C * m_x)).eval();
        if(!next_x.allFinite())
            return ctrlpp::unexpected(
                luenberger_update_error::non_finite_result);
        m_x = std::move(next_x);
        return {};
    }

    auto state() const -> const state_vector_t& { return m_x; }

    /// @brief Report whether the carried estimate is still degraded from an
    /// earlier step. Latches until `reset`; a rejected measurement does not set
    /// it.
    auto health() const -> luenberger_health { return m_health; }

    void set_gain(const gain_matrix_t& L) { m_L = L; }
    void set_model(system_t sys) { m_sys = std::move(sys); }

    /// @brief Seed a fresh estimate. This is the one operation that clears a
    /// latched degraded status, because it replaces the very state the status
    /// describes.
    void reset(const state_vector_t& x0)
    {
        m_x = x0;
        m_health = luenberger_health::ok;
    }

private:
    /// @brief Classify a step's operands without touching a single member.
    ///
    /// The order is the severity order documented on
    /// `luenberger_update_error`: the carried estimate first, the supplied
    /// measurement last. The cost is one finiteness scan of each operand --
    /// NX + NY reads, no branches on data and no allocation, both dimensions
    /// being compile-time constants.
    auto check_step(const output_vector_t& z) const -> ctrlpp::expected<void, luenberger_update_error>
    {
        if(!m_x.allFinite())
            return ctrlpp::unexpected(luenberger_update_error::non_finite_state);
        if(!z.allFinite())
            return ctrlpp::unexpected(luenberger_update_error::non_finite_measurement);
        return {};
    }

    /// @brief Latch the persistent status for the fault that describes the
    /// carried estimate, and pass the fault through so the caller still receives
    /// the specific cause on the failure channel.
    auto latch_health(luenberger_update_error fault) -> luenberger_update_error
    {
        if(fault != luenberger_update_error::non_finite_measurement)
            m_health = luenberger_health::non_finite_estimate;
        return fault;
    }

    system_t m_sys;
    gain_matrix_t m_L;
    state_vector_t m_x;
    luenberger_health m_health{luenberger_health::ok};
};

static_assert(ObserverPolicy<luenberger_observer<double, 2, 1, 1>>);
static_assert(!CovarianceObserver<luenberger_observer<double, 2, 1, 1>>);

}

#endif
