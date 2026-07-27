#ifndef HPP_GUARD_CTRLPP_NMPC_H
#define HPP_GUARD_CTRLPP_NMPC_H

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/nmpc_config.h"
#include "ctrlpp/mpc/diagnostics.h"
#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/mpc/nlp_formulation.h"

#include <Eigen/Dense>

#include <span>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>
#include <cstddef>
#include <utility>
#include <optional>
#include <algorithm>

namespace ctrlpp
{

/// @brief Runtime-horizon nonlinear MPC — the opt-in soft-RT / non-RT path.
///
/// The horizon is a runtime field (`nmpc_config::horizon`), so the decision
/// dimension is only known at construction and the solve binds a runtime-erased
/// `nlp_problem<Scalar>`, forcing the solver's allocating dynamic path. This is
/// the flexible escape hatch, NOT the default: the public `nmpc` name resolves to
/// the compile-time-horizon `nmpc_static` (allocation-free at argmin's fixed-N
/// floor). Select `nmpc_dynamic` explicitly when the horizon must vary at runtime
/// and soft-RT latency is acceptable. See the RT-safety matrix.
template <typename Scalar, std::size_t NX, std::size_t NU, nlp_solver Solver, dynamics_model<Scalar, NX, NU> Dynamics, std::size_t NC = 0, std::size_t NTC = 0>
class nmpc_dynamic
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);
    static constexpr int nc = static_cast<int>(NC);
    static constexpr int ntc = static_cast<int>(NTC);

public:
    /// @brief Validating factory; the only construction path on an
    /// exception-free build. Builds the solver with its own defaults and chains
    /// into the solver-taking overload below.
    ///
    /// Rejections, checked in order:
    ///  * horizon <= 0                   -> controller_construction_error::non_positive_horizon
    ///  * horizon above the representable
    ///    bound of the derived dimensions -> controller_construction_error::horizon_overflow
    [[nodiscard]] static auto create(Dynamics dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config)
        -> expected<nmpc_dynamic, controller_construction_error>
    {
        return create(std::move(dynamics), config, Solver{});
    }

    /// @brief Validating factory taking a caller-supplied, pre-configured
    /// solver. Both the dynamics and the solver are moved in before the NLP is
    /// posed, so the caller's solver settings govern setup.
    ///
    /// The horizon is the only runtime quantity that scales the posed problem,
    /// and it is validated here, before any dimension product is formed and
    /// before any storage is reserved. Rejections, checked in order:
    ///  * horizon <= 0             -> controller_construction_error::non_positive_horizon
    ///  * horizon > horizon_bound   -> controller_construction_error::horizon_overflow
    ///
    /// The overflow bound is a representability condition on the horizon's own
    /// type, not a chosen ceiling. At the worst-case configuration (path slack,
    /// terminal slack and rate bounds all present) the horizon N scales two
    /// dimensions:
    ///   decision vector : (N+1)*nx + N*nu + N*nc + ntc
    ///                                       = N*(nx + nu + nc) + nx + ntc
    ///   constraint rows : (N+1)*nx + 2*N*nu + N*nc + ntc
    ///                                       = N*(nx + 2*nu + nc) + nx + ntc
    /// so the largest per-step contribution is nx + 2*nu + nc and the
    /// horizon-independent part is nx + ntc. Both products therefore stay
    /// representable exactly when
    ///   horizon <= (max<int> - (nx + ntc)) / (nx + 2*nu + nc).
    ///
    /// Two alternatives are deliberately not implemented. Validating at the
    /// first solve would surface a configuration error at the first control
    /// step, the worst possible moment. Clamping the horizon to one would turn a
    /// caller mistake into a silently different controller.
    [[nodiscard]] static auto create(Dynamics dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config, Solver solver)
        -> expected<nmpc_dynamic, controller_construction_error>
    {
        if(config.horizon <= 0)
            return unexpected(controller_construction_error::non_positive_horizon);
        if(config.horizon > horizon_bound())
            return unexpected(controller_construction_error::horizon_overflow);

        return nmpc_dynamic{unchecked_t{}, std::move(dynamics), config, std::move(solver)};
    }

private:
    /// @brief Tag selecting the non-validating constructor reserved for `create`.
    struct unchecked_t
    {
        explicit unchecked_t() = default;
    };

    /// @brief Largest horizon whose derived decision and constraint dimensions
    /// are still representable in the horizon's own type. See `create` for
    /// the derivation; this forms no product of its own.
    [[nodiscard]] static auto horizon_bound() -> int
    {
        constexpr int per_step = nx + 2 * nu + nc;
        constexpr int constant_dimensions = nx + ntc;
        return (std::numeric_limits<int>::max() - constant_dimensions) / per_step;
    }

    /// @brief Construct from a configuration already validated by `create`.
    nmpc_dynamic(unchecked_t, Dynamics dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config, Solver solver)
        : m_dynamics{std::move(dynamics)}
        , m_config{config}
        , m_N{config.horizon}
        , m_has_path_slack{config.soft_constraints && NC > 0 && config.path_constraint.has_value()}
        , m_has_term_slack{config.soft_constraints && NTC > 0 && config.terminal_constraint.has_value()}
        , m_num_path_slack{m_has_path_slack ? m_N * nc : 0}
        , m_num_term_slack{m_has_term_slack ? ntc : 0}
        , m_num_vars{(m_N + 1) * nx + m_N * nu + m_num_path_slack + m_num_term_slack}
        , m_state{std::make_shared<nmpc_formulation_state<Scalar, NX, NU>>()}
        , m_solver{std::move(solver)}
    {
        m_state->x_ref.resize(static_cast<std::size_t>(m_N + 1), Vector<Scalar, NX>::Zero());
        m_problem = std::make_unique<nlp_problem<Scalar>>(
            detail::build_nmpc_problem<Scalar, NX, NU, NC, NTC>(m_dynamics, m_config, m_state));
        m_setup_failed = !detail::setup_nlp_solver(m_solver, *m_problem);
        m_warm_z = Eigen::VectorX<Scalar>::Zero(m_num_vars);
    }

public:
    // Move is correct-by-default: `m_problem` lives behind a `unique_ptr` (a
    // stable heap address), so moving relocates only the owning pointer while the
    // pointed-to problem stays put. The solver's bridge caches `&m_problem` by
    // reference (argmin_problem.h:24-27), so that cached pointer stays valid
    // across a defaulted move. `m_state` is already a heap-stable shared_ptr.
    nmpc_dynamic(nmpc_dynamic&&) = default;
    nmpc_dynamic& operator=(nmpc_dynamic&&) = default;
    ~nmpc_dynamic() = default;

    // Copy is an independent fork (the solver's copy is copy-deleted-by-argmin
    // and re-emplaced lazily, so aliasing would be unsafe anyway): deep-copy the
    // formulation state, rebuild the problem bound to that fresh state, then
    // re-setup the forked solver against the new problem. The two controllers
    // share no mutable state and solve independently.
    nmpc_dynamic(const nmpc_dynamic& other)
        : m_dynamics{other.m_dynamics}
        , m_config{other.m_config}
        , m_N{other.m_N}
        , m_has_path_slack{other.m_has_path_slack}
        , m_has_term_slack{other.m_has_term_slack}
        , m_num_path_slack{other.m_num_path_slack}
        , m_num_term_slack{other.m_num_term_slack}
        , m_num_vars{other.m_num_vars}
        , m_state{std::make_shared<nmpc_formulation_state<Scalar, NX, NU>>(*other.m_state)}
        , m_problem{std::make_unique<nlp_problem<Scalar>>(
              detail::build_nmpc_problem<Scalar, NX, NU, NC, NTC>(m_dynamics, m_config, m_state))}
        , m_solver{other.m_solver}
        , m_warm_z{other.m_warm_z}
        , m_last_solution{other.m_last_solution}
        , m_last_diagnostics{other.m_last_diagnostics}
        , m_u_prev{other.m_u_prev}
        , m_has_solution{other.m_has_solution}
    {
        m_setup_failed = !detail::setup_nlp_solver(m_solver, *m_problem);
    }

    nmpc_dynamic& operator=(const nmpc_dynamic& other)
    {
        nmpc_dynamic tmp{other};
        *this = std::move(tmp);
        return *this;
    }

    // Unified soft-constraint / failure contract, identical to mpc (see mpc.h).
    //
    // solve() returns ctrlpp::expected<solve_output<Scalar, NU>, solver_error>.
    //   * SUCCESS branch: a usable control input reached through `->input`, plus a
    //     soft solve_result_status (optimal->converged,
    //     solved_inaccurate->solved_inaccurate,
    //     max_iterations/time_limit->budget_exhausted). No implicit conversion to
    //     Vector, so the status is never silently dropped.
    //   * ERROR branch: a hard failure (infeasible/invalid_problem/setup_incomplete/
    //     invalid_backend_result) with NO input. On error the internal m_u_prev is
    //     NOT updated and no hidden fallback input is applied. Use
    //     set_applied_input to record the input the caller actually commanded.
    expected<solve_output<Scalar, NU>, solver_error> solve(const Vector<Scalar, NX>& x0)
    {
        for(auto& ref : m_state->x_ref)
            ref.setZero();
        return solve_impl(x0);
    }

    expected<solve_output<Scalar, NU>, solver_error> solve(const Vector<Scalar, NX>& x0, const Vector<Scalar, NX>& x_ref)
    {
        for(auto& ref : m_state->x_ref)
            ref = x_ref;
        return solve_impl(x0);
    }

    expected<solve_output<Scalar, NU>, solver_error> solve(const Vector<Scalar, NX>& x0, std::span<const Vector<Scalar, NX>> x_ref)
    {
        // An empty reference span has no value to back-fill from; reject it via
        // the error branch rather than solving against stale references.
        if(x_ref.empty())
            return unexpected(solver_error::invalid_problem);
        const auto len = std::min(x_ref.size(), m_state->x_ref.size());
        for(std::size_t k = 0; k < len; ++k)
            m_state->x_ref[k] = x_ref[k];
        for(std::size_t k = len; k < m_state->x_ref.size(); ++k)
            m_state->x_ref[k] = x_ref.back();
        return solve_impl(x0);
    }

    /// @brief Record the control input the caller actually commanded.
    ///
    /// The internal m_u_prev is updated only on a successful solve; this accessor
    /// lets the caller keep it consistent when it commands a different input.
    void set_applied_input(const Vector<Scalar, NU>& u) { m_u_prev = u; }

    // Guarded: returns the error branch before the first valid solve, so a caller
    // can never read stale or default-initialized solution data.
    expected<std::pair<std::vector<Vector<Scalar, NX>>, std::vector<Vector<Scalar, NU>>>, solver_error> trajectory() const
    {
        if(!m_has_solution)
            return unexpected(solver_error::setup_incomplete);

        // The slices below need no length check of their own: m_has_solution is
        // set only by finish_solve, which stores the primal only after checking
        // that it covers the problem dimension these offsets are derived from.
        std::vector<Vector<Scalar, NX>> states;
        std::vector<Vector<Scalar, NU>> inputs;
        states.reserve(static_cast<std::size_t>(m_N + 1));
        inputs.reserve(static_cast<std::size_t>(m_N));

        for(int k = 0; k <= m_N; ++k)
            states.push_back(m_last_solution.segment(k * nx, nx));
        const int u_offset = (m_N + 1) * nx;
        for(int k = 0; k < m_N; ++k)
            inputs.push_back(m_last_solution.segment(u_offset + k * nu, nu));

        return std::pair{std::move(states), std::move(inputs)};
    }

    mpc_diagnostics<Scalar> diagnostics() const { return m_last_diagnostics; }
    const nlp_problem<Scalar>& problem() const { return *m_problem; }
    const Eigen::VectorX<Scalar>& last_solution() const { return m_last_solution; }

private:
    expected<solve_output<Scalar, NU>, solver_error> solve_impl(const Vector<Scalar, NX>& x0)
    {
        if(m_setup_failed)
            return unexpected(solver_error::setup_incomplete);

        m_state->x0 = x0;
        m_state->u_prev = m_u_prev;

        nlp_update<Scalar> update;
        update.x0 = m_warm_z;

        auto result = m_solver.solve(update);
        populate_diagnostics(result);

        // WIDENED accept-set: budget-limited iterates (max_iterations/time_limit)
        // now reach the caller on the SUCCESS branch tagged budget_exhausted.
        switch(result.status)
        {
        case solve_status::optimal:
            return finish_solve(result.x, solve_result_status::converged);
        case solve_status::solved_inaccurate:
            return finish_solve(result.x, solve_result_status::solved_inaccurate);
        case solve_status::max_iterations:
        case solve_status::time_limit:
            return finish_solve(result.x, solve_result_status::budget_exhausted);
        case solve_status::infeasible:
            return unexpected(solver_error::infeasible);
        case solve_status::unbounded:
        case solve_status::non_convex:
        case solve_status::error:
        default:
            return unexpected(solver_error::invalid_problem);
        }
    }

    /// @brief Consume an accepted backend result: validate its reported shape,
    /// then store it, diagnose it, shift the warm start and extract the input.
    ///
    /// The accept-set above is decided purely from the status the backend
    /// reports, and a status is not a shape. Every read below -- the constraint
    /// diagnostics' raw pointer maps, the warm-start shift, and the first-input
    /// slice -- sits at an offset derived from the horizon, so all of them are
    /// covered by this one comparison against the problem dimension the
    /// controller derived at construction from its horizon and its state, input,
    /// path-slack and terminal-slack contributions. Nothing is stored and no
    /// member is touched when it fails.
    [[nodiscard]] auto finish_solve(const Eigen::VectorX<Scalar>& z, solve_result_status status) -> expected<solve_output<Scalar, NU>, solver_error>
    {
        if(z.size() < static_cast<Eigen::Index>(m_num_vars))
        {
            // Correct the diagnostics the backend's own status just populated, so
            // a reader there is not told the solve was optimal when its answer
            // was discarded.
            m_last_diagnostics.status = solve_status::invalid_backend_result;
            return unexpected(solver_error::invalid_backend_result);
        }

        m_last_solution = z;
        m_has_solution = true;
        populate_constraint_diagnostics(z);
        shift_warm_start(z);
        return solve_output<Scalar, NU>{.input = extract_first_input(z), .status = status};
    }

    void populate_diagnostics(const nlp_result<Scalar>& result)
    {
        m_last_diagnostics = mpc_diagnostics<Scalar>{.status = result.status,
                                                     .iterations = result.iterations,
                                                     .solve_time = result.solve_time,
                                                     .cost = result.objective,
                                                     .primal_residual = result.primal_residual,
                                                     .dual_residual = Scalar{0},
                                                     .max_constraint_violation = result.primal_residual,
                                                     .max_path_constraint_violation = Scalar{0},
                                                     .max_terminal_constraint_violation = Scalar{0},
                                                     .total_slack = Scalar{0}};
    }

    void populate_constraint_diagnostics(const Eigen::VectorX<Scalar>& z)
    {
        const int u_offset = (m_N + 1) * nx;
        compute_path_constraint_violation(z, u_offset);
        compute_terminal_constraint_violation(z);
        compute_total_slack(z, u_offset);
    }

    void compute_path_constraint_violation(const Eigen::VectorX<Scalar>& z, int u_offset)
    {
        if constexpr(NC > 0)
        {
            if(m_config.path_constraint)
            {
                const auto& g = *m_config.path_constraint;
                Scalar max_viol{0};

                for(int k = 0; k < m_N; ++k)
                {
                    Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + k * nx);
                    Eigen::Map<const Vector<Scalar, NU>> uk(z.data() + u_offset + k * nu);
                    Vector<Scalar, NC> gk = g(xk, uk);

                    for(int i = 0; i < nc; ++i)
                        max_viol = std::max(max_viol, gk[i]);
                }

                m_last_diagnostics.max_path_constraint_violation = max_viol;
            }
        }
    }

    void compute_terminal_constraint_violation(const Eigen::VectorX<Scalar>& z)
    {
        if constexpr(NTC > 0)
        {
            if(m_config.terminal_constraint)
            {
                const auto& h = *m_config.terminal_constraint;
                Eigen::Map<const Vector<Scalar, NX>> xN(z.data() + m_N * nx);
                Vector<Scalar, NTC> hN = h(xN);

                Scalar max_viol{0};
                for(int i = 0; i < ntc; ++i)
                    max_viol = std::max(max_viol, hN[i]);

                m_last_diagnostics.max_terminal_constraint_violation = max_viol;
            }
        }
    }

    void compute_total_slack(const Eigen::VectorX<Scalar>& z, int u_offset)
    {
        Scalar total{0};
        if(m_has_path_slack)
        {
            const int slack_off = u_offset + m_N * nu;
            for(int k = 0; k < m_N; ++k)
            {
                Eigen::Map<const Vector<Scalar, NC>> sk(z.data() + slack_off + k * nc);
                total += sk.sum();
            }
        }
        if(m_has_term_slack)
        {
            const int term_off = u_offset + m_N * nu + m_num_path_slack;
            Eigen::Map<const Vector<Scalar, NTC>> st(z.data() + term_off);
            total += st.sum();
        }
        m_last_diagnostics.total_slack = total;
    }

    void shift_warm_start(const Eigen::VectorX<Scalar>& sol)
    {
        m_warm_z = sol;
        shift_warm_states(sol);
        shift_warm_inputs(sol);
        shift_warm_path_slack(sol);
        zero_warm_terminal_slack();
    }

    void shift_warm_states(const Eigen::VectorX<Scalar>& sol)
    {
        for(int k = 0; k < m_N; ++k)
            m_warm_z.segment(k * nx, nx) = sol.segment((k + 1) * nx, nx);
    }

    void shift_warm_inputs(const Eigen::VectorX<Scalar>& sol)
    {
        const int u_offset = (m_N + 1) * nx;
        for(int k = 0; k < m_N - 1; ++k)
            m_warm_z.segment(u_offset + k * nu, nu) = sol.segment(u_offset + (k + 1) * nu, nu);
    }

    void shift_warm_path_slack(const Eigen::VectorX<Scalar>& sol)
    {
        if(!m_has_path_slack)
            return;
        const int u_offset = (m_N + 1) * nx;
        const int slack_off = u_offset + m_N * nu;
        for(int k = 0; k < m_N - 1; ++k)
            m_warm_z.segment(slack_off + k * nc, nc) = sol.segment(slack_off + (k + 1) * nc, nc);
        m_warm_z.segment(slack_off + (m_N - 1) * nc, nc).setZero();
    }

    void zero_warm_terminal_slack()
    {
        if(!m_has_term_slack)
            return;
        const int u_offset = (m_N + 1) * nx;
        const int term_off = u_offset + m_N * nu + m_num_path_slack;
        m_warm_z.segment(term_off, ntc).setZero();
    }

    auto extract_first_input(const Eigen::VectorX<Scalar>& sol) -> Vector<Scalar, NU>
    {
        const int u_offset = (m_N + 1) * nx;
        Vector<Scalar, NU> u0 = sol.segment(u_offset, nu);
        m_u_prev = u0;
        return u0;
    }

    Dynamics m_dynamics;
    nmpc_config<Scalar, NX, NU, NC, NTC> m_config;
    int m_N;
    bool m_has_path_slack;
    bool m_has_term_slack;
    int m_num_path_slack;
    int m_num_term_slack;
    int m_num_vars;

    std::shared_ptr<nmpc_formulation_state<Scalar, NX, NU>> m_state;
    // Held behind a stable heap address so a defaulted move does not relocate the
    // object the solver's bridge points at by reference (argmin_problem.h:24-27).
    std::unique_ptr<nlp_problem<Scalar>> m_problem;
    Solver m_solver{};

    Eigen::VectorX<Scalar> m_warm_z;
    Eigen::VectorX<Scalar> m_last_solution;
    mpc_diagnostics<Scalar> m_last_diagnostics{};
    Vector<Scalar, NU> m_u_prev{Vector<Scalar, NU>::Zero()};
    bool m_setup_failed{false};
    bool m_has_solution{false};
};

/// @brief Compile-time-horizon nonlinear MPC (additive Route A / SEED-002).
///
/// Parallel to the runtime-horizon `nmpc_dynamic`, but the horizon NH is a template
/// parameter, so the decision dimension NV = (NH+1)*NX + NH*NU is a compile-time
/// constant threaded through `build_nmpc_problem_static` into an
/// `nlp_problem_static<Scalar, NV>` and a compile-time-N `argmin_solver`. On
/// argmin's fixed-N NW-SQP floor the steady-state solve is allocation-free with
/// stock Eigen; `nmpc_static_nomalloc_test` is the proof.
///
/// This is the slack-free (hard-constraint) cut: config.horizon must equal the
/// compile-time NH and the formulation must not introduce slack decision
/// variables, so n_vars == NV exactly (enforced by build_nmpc_problem_static).
///
/// Solver-generic by design: the caller supplies the static solver type already
/// parameterized on the matching NV (e.g. `argmin_solver<Scalar, argmin_nw_sqp,
/// true, NV>`), so this header pulls in no argmin dependency and the argmin-off
/// build compiles unchanged. The public `nmpc` name is an alias to `nmpc_static`
/// (defined below the class), so the compile-time-horizon, allocation-free path is
/// the DEFAULT; the runtime-horizon class is the opt-in `nmpc_dynamic`.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NH, typename Solver, dynamics_model<Scalar, NX, NU> Dynamics, std::size_t NC = 0, std::size_t NTC = 0>
class nmpc_static
{
    static_assert(NH > 0, "Horizon NH must be positive: a zero horizon leaves no input to apply and makes the horizon - 1 warm-start shift and the input offset into the primal read outside the decision vector");

    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);

public:
    static constexpr int horizon = static_cast<int>(NH);
    static constexpr int problem_dimension = static_cast<int>((NH + 1) * NX + NH * NU);

    using problem_type = nlp_problem_static<Scalar, problem_dimension>;

    nmpc_static(Dynamics dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config)
        : nmpc_static{std::move(dynamics), config, Solver{}}
    {}

    nmpc_static(Dynamics dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config, Solver solver)
        : m_dynamics{std::move(dynamics)}
        , m_config{config}
        , m_state{std::make_shared<nmpc_formulation_state<Scalar, NX, NU>>()}
        , m_solver{std::move(solver)}
    {
        m_state->x_ref.resize(static_cast<std::size_t>(horizon) + 1, Vector<Scalar, NX>::Zero());
        build_problem_and_setup();
        // Pre-size the reused solve buffers ONCE (construction may allocate), so
        // the steady-state assignments below hit the same-size fast path and the
        // hot solve loop stays allocation-free.
        m_warm_z = Eigen::VectorX<Scalar>::Zero(problem_dimension);
        m_update.x0 = Eigen::VectorX<Scalar>::Zero(problem_dimension);
        m_last_solution = Eigen::VectorX<Scalar>::Zero(problem_dimension);
    }

    // Move is correct-by-default for the same reason as runtime-horizon nmpc:
    // `m_problem` lives behind a unique_ptr (stable heap address), so the solver
    // bridge's by-reference back-pointer stays valid across a defaulted move.
    nmpc_static(nmpc_static&&) = default;
    nmpc_static& operator=(nmpc_static&&) = default;
    ~nmpc_static() = default;

    // Copy is an independent fork: deep-copy the formulation state, rebuild the
    // problem bound to that fresh state, then re-setup the forked solver.
    nmpc_static(const nmpc_static& other)
        : m_dynamics{other.m_dynamics}
        , m_config{other.m_config}
        , m_state{std::make_shared<nmpc_formulation_state<Scalar, NX, NU>>(*other.m_state)}
        , m_solver{other.m_solver}
        , m_warm_z{other.m_warm_z}
        , m_last_solution{other.m_last_solution}
        , m_update{other.m_update}
        , m_last_diagnostics{other.m_last_diagnostics}
        , m_u_prev{other.m_u_prev}
        , m_has_solution{other.m_has_solution}
    {
        build_problem_and_setup();
    }

    nmpc_static& operator=(const nmpc_static& other)
    {
        nmpc_static tmp{other};
        *this = std::move(tmp);
        return *this;
    }

    expected<solve_output<Scalar, NU>, solver_error> solve(const Vector<Scalar, NX>& x0)
    {
        for(auto& ref : m_state->x_ref)
            ref.setZero();
        return solve_impl(x0);
    }

    expected<solve_output<Scalar, NU>, solver_error> solve(const Vector<Scalar, NX>& x0, const Vector<Scalar, NX>& x_ref)
    {
        for(auto& ref : m_state->x_ref)
            ref = x_ref;
        return solve_impl(x0);
    }

    void set_applied_input(const Vector<Scalar, NU>& u) { m_u_prev = u; }

    mpc_diagnostics<Scalar> diagnostics() const { return m_last_diagnostics; }
    const problem_type& problem() const { return *m_problem; }
    const Eigen::VectorX<Scalar>& last_solution() const { return m_last_solution; }

private:
    // Pose the compile-time-dimension problem and set the solver up against it.
    //
    // The formulation factory rejects a configuration whose runtime horizon
    // disagrees with the compile-time NH, or that would introduce slack decision
    // variables, and it does so before building any of the dependent problem
    // data. That rejection is unconditional (a release build no longer steps
    // over it), and it is reported to the caller through the controller's
    // existing error channel: the setup flag, which turns every subsequent
    // solve() into solver_error::setup_incomplete. On rejection the problem
    // holder is still populated with an empty problem so problem() has something
    // valid to reference.
    void build_problem_and_setup()
    {
        auto problem = detail::build_nmpc_problem_static<Scalar, NX, NU, NH, NC, NTC>(m_dynamics, m_config, m_state);
        if(!problem.has_value())
        {
            m_problem = std::make_unique<problem_type>();
            m_setup_failed = true;
            return;
        }

        m_problem = std::make_unique<problem_type>(*std::move(problem));
        m_setup_failed = !setup_solver();
    }

    // Solver setup without the runtime-erased setup_nlp_solver helper (that helper
    // is typed on nlp_problem<Scalar>; the static path binds nlp_problem_static).
    // Kept solver-generic: fallible try_setup when available, classic setup
    // otherwise.
    [[nodiscard]] bool setup_solver()
    {
        if constexpr(requires { m_solver.try_setup(*m_problem); })
            return m_solver.try_setup(*m_problem).has_value();
        else
        {
            m_solver.setup(*m_problem);
            return true;
        }
    }

    expected<solve_output<Scalar, NU>, solver_error> solve_impl(const Vector<Scalar, NX>& x0)
    {
        if(m_setup_failed)
            return unexpected(solver_error::setup_incomplete);

        m_state->x0 = x0;
        m_state->u_prev = m_u_prev;

        // Same-size assignment into the pre-sized member: no reallocation once the
        // constructor has sized x0 to the compile-time NV.
        m_update.x0 = m_warm_z;

        const nlp_result<Scalar> result = dispatch_solve();
        populate_diagnostics(result);

        switch(result.status)
        {
        case solve_status::optimal:
            return finish_solve(solve_result_status::converged);
        case solve_status::solved_inaccurate:
            return finish_solve(solve_result_status::solved_inaccurate);
        case solve_status::max_iterations:
        case solve_status::time_limit:
            return finish_solve(solve_result_status::budget_exhausted);
        case solve_status::infeasible:
            return unexpected(solver_error::infeasible);
        case solve_status::unbounded:
        case solve_status::non_convex:
        case solve_status::error:
        default:
            return unexpected(solver_error::invalid_problem);
        }
    }

    // Detection: the static (strict-zero) argmin bridge exposes a zero-alloc
    // solve_into; the runtime-erased solvers expose only solve(). Prefer
    // solve_into so the steady-state hot path allocates no dynamic decision
    // vector (solve() must return nlp_result::x by value, which heap-allocates).
    static constexpr bool solver_has_solve_into =
        requires(Solver& s, const nlp_update<Scalar>& u, Eigen::VectorX<Scalar>& out) { s.solve_into(u, out); };

    // Run the solver, always leaving the primal in the pre-sized m_last_solution
    // and returning the diagnostics. The solve_into branch writes m_last_solution
    // in place (no allocation); the solve() fallback copies its by-value primal.
    nlp_result<Scalar> dispatch_solve()
    {
        if constexpr(solver_has_solve_into)
            return m_solver.solve_into(m_update, m_last_solution);
        else
        {
            auto result = m_solver.solve(m_update);
            m_last_solution = result.x;
            return result;
        }
    }

    /// @brief Consume an accepted backend result: validate the shape of the
    /// primal the dispatch left behind, then shift the warm start and extract the
    /// input.
    ///
    /// One comparison covers both dispatch shapes. The write-into branch fills a
    /// buffer this class pre-sized to the compile-time dimension, so there the
    /// check is near-always trivially true; the fallback branch assigns a
    /// by-value primal of whatever length the backend chose, which resizes the
    /// buffer and is exactly the case that must not reach the slices below.
    /// Checking m_last_solution after the dispatch, rather than the result of
    /// either branch, is what makes the single check sufficient.
    [[nodiscard]] auto finish_solve(solve_result_status status) -> expected<solve_output<Scalar, NU>, solver_error>
    {
        if(m_last_solution.size() < static_cast<Eigen::Index>(problem_dimension))
        {
            // Correct the diagnostics the backend's own status just populated, so
            // a reader there is not told the solve was optimal when its answer
            // was discarded.
            m_last_diagnostics.status = solve_status::invalid_backend_result;
            return unexpected(solver_error::invalid_backend_result);
        }

        m_has_solution = true;
        shift_warm_start(m_last_solution);
        return solve_output<Scalar, NU>{.input = extract_first_input(m_last_solution), .status = status};
    }

    void populate_diagnostics(const nlp_result<Scalar>& result)
    {
        m_last_diagnostics = mpc_diagnostics<Scalar>{.status = result.status,
                                                     .iterations = result.iterations,
                                                     .solve_time = result.solve_time,
                                                     .cost = result.objective,
                                                     .primal_residual = result.primal_residual,
                                                     .dual_residual = Scalar{0},
                                                     .max_constraint_violation = result.primal_residual,
                                                     .max_path_constraint_violation = Scalar{0},
                                                     .max_terminal_constraint_violation = Scalar{0},
                                                     .total_slack = Scalar{0}};
    }

    // Warm-start shift for the slack-free formulation: advance the state and input
    // blocks by one node; the tail input is held. In-place segment writes on the
    // pre-sized m_warm_z, so no allocation.
    void shift_warm_start(const Eigen::VectorX<Scalar>& sol)
    {
        m_warm_z = sol;
        for(int k = 0; k < horizon; ++k)
            m_warm_z.segment(k * nx, nx) = sol.segment((k + 1) * nx, nx);
        const int u_offset = (horizon + 1) * nx;
        for(int k = 0; k < horizon - 1; ++k)
            m_warm_z.segment(u_offset + k * nu, nu) = sol.segment(u_offset + (k + 1) * nu, nu);
    }

    auto extract_first_input(const Eigen::VectorX<Scalar>& sol) -> Vector<Scalar, NU>
    {
        const int u_offset = (horizon + 1) * nx;
        Vector<Scalar, NU> u0 = sol.segment(u_offset, nu);
        m_u_prev = u0;
        return u0;
    }

    Dynamics m_dynamics;
    nmpc_config<Scalar, NX, NU, NC, NTC> m_config;

    std::shared_ptr<nmpc_formulation_state<Scalar, NX, NU>> m_state;
    // Held behind a stable heap address so a defaulted move does not relocate the
    // object the solver's bridge points at by reference.
    std::unique_ptr<problem_type> m_problem;
    Solver m_solver{};

    Eigen::VectorX<Scalar> m_warm_z;
    Eigen::VectorX<Scalar> m_last_solution;
    nlp_update<Scalar> m_update{};
    mpc_diagnostics<Scalar> m_last_diagnostics{};
    Vector<Scalar, NU> m_u_prev{Vector<Scalar, NU>::Zero()};
    bool m_setup_failed{false};
    bool m_has_solution{false};
};

/// @brief Public default NMPC name: the compile-time-horizon, allocation-free
/// `nmpc_static`. Pinning the horizon NH at compile time makes the steady-state
/// solve allocation-free at argmin's fixed-N floor (see `nmpc_static`); use
/// `nmpc_dynamic` when the horizon must vary at runtime (soft-RT). This re-points
/// the public `nmpc` name from the former runtime-horizon class (now
/// `nmpc_dynamic`) to the static default — a deliberate pre-1.0 public-API-shape
/// change.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NH, typename Solver, dynamics_model<Scalar, NX, NU> Dynamics, std::size_t NC = 0, std::size_t NTC = 0>
using nmpc = nmpc_static<Scalar, NX, NU, NH, Solver, Dynamics, NC, NTC>;

}

#endif
