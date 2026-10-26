#include "ctrlpp/mpc/argmin_qp_solver.h"

#include <catch2/catch_test_macros.hpp>

// Catch2 exits non-zero from a binary with no registered test, so this case stays outside the gate.
TEST_CASE("mpc and mhe on argmin_qp_solver gate", "[mpc][mhe][argmin][qp]")
{
#if defined(CTRLPP_HAS_ARGMIN_QP)
    SUCCEED("argmin QP header present: the mpc and mhe cases below run");
#else
    SUCCEED("argmin QP header absent: the mpc and mhe cases are compiled out");
#endif
}

#if defined(CTRLPP_HAS_ARGMIN_QP)

#include "hardening_helpers.h"

#include "ctrlpp/mhe.h"
#include "ctrlpp/mpc.h"

#include "ctrlpp/control/dare.h"

#include "ctrlpp/estimation/kalman.h"

#include "ctrlpp/mpc/qp_formulation.h"

#if defined(CTRLPP_HAS_OSQP)
#include "ctrlpp/mpc/osqp_solver.h"
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cmath>
#include <limits>
#include <random>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <optional>
#include <algorithm>
#include <type_traits>

static_assert(std::is_move_constructible_v<ctrlpp::argmin_qp_solver>, "mpc and mhe hold their solver by value");

namespace
{

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr std::size_t NY = 1;
constexpr std::size_t N = 5;
constexpr double dt = 0.1;

auto make_double_integrator() -> ctrlpp::discrete_state_space<double, NX, NU, NX>
{
    Eigen::Matrix2d A;
    A << 1.0, dt, 0.0, 1.0;
    Eigen::Vector2d B;
    B << 0.5 * dt * dt, dt;
    Eigen::Matrix2d C = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 2, 1> D = Eigen::Matrix<double, 2, 1>::Zero();
    return {A, B, C, D};
}

struct linear_dynamics
{
    auto operator()(const ctrlpp::Vector<double, NX>& x, const ctrlpp::Vector<double, NU>& u) const -> ctrlpp::Vector<double, NX>
    {
        ctrlpp::Vector<double, NX> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    }

    auto jacobian_x(const ctrlpp::Vector<double, NX>& /*x*/, const ctrlpp::Vector<double, NU>& /*u*/) const -> ctrlpp::Matrix<double, NX, NX>
    {
        ctrlpp::Matrix<double, NX, NX> F;
        F << 1.0, dt, 0.0, 1.0;
        return F;
    }

    auto jacobian_u(const ctrlpp::Vector<double, NX>& /*x*/, const ctrlpp::Vector<double, NU>& /*u*/) const -> ctrlpp::Matrix<double, NX, NU>
    {
        ctrlpp::Matrix<double, NX, NU> G;
        G << 0.5 * dt * dt, dt;
        return G;
    }
};

struct position_measurement
{
    auto operator()(const ctrlpp::Vector<double, NX>& x) const -> ctrlpp::Vector<double, NY>
    {
        return (ctrlpp::Vector<double, NY>() << x(0)).finished();
    }

    auto jacobian(const ctrlpp::Vector<double, NX>& /*x*/) const -> ctrlpp::Matrix<double, NY, NX>
    {
        return (ctrlpp::Matrix<double, NY, NX>() << 1.0, 0.0).finished();
    }
};

using argmin_mpc = ctrlpp::mpc<double, NX, NU, ctrlpp::argmin_qp_solver>;
using argmin_mhe = ctrlpp::mhe<double, NX, NU, NY, N, ctrlpp::argmin_qp_solver, linear_dynamics, position_measurement>;

}

TEST_CASE("mpc on argmin_qp_solver regulates the double integrator", "[mpc][argmin][qp]")
{
    const auto sys = make_double_integrator();
    const ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 10,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
    };

    auto controller = ctrlpp::test::constructed(argmin_mpc::create(sys, cfg, ctrlpp::argmin_qp_solver{ctrlpp::qp_preset::speed}));

    Eigen::Vector2d x{1.0, 0.0};
    for(std::size_t step = 0; step < 50; ++step)
    {
        const auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = sys.A * x + sys.B * u->input;
    }
    CHECK(x.norm() < 0.1);
}

TEST_CASE("mhe on argmin_qp_solver solves its window after warm-up", "[mhe][argmin][qp]")
{
    ctrlpp::mhe_config<double, NX, NU, NY, N> cfg;
    cfg.Q = ctrlpp::Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = ctrlpp::Matrix<double, NY, NY>::Identity() * 0.1;
    cfg.P0 = ctrlpp::Matrix<double, NX, NX>::Identity() * 10.0;

    auto estimator = ctrlpp::test::constructed(argmin_mhe::create(linear_dynamics{}, position_measurement{}, cfg));

    ctrlpp::Vector<double, NX> x_true{0.0, 1.0};
    const ctrlpp::Vector<double, NU> u = ctrlpp::Vector<double, NU>::Zero();
    for(std::size_t i = 0; i < N + 10; ++i)
    {
        x_true(0) += dt * x_true(1);
        estimator.predict(u);
        REQUIRE(estimator.update((ctrlpp::Vector<double, NY>() << x_true(0)).finished()).has_value());

        // The window QP first runs on the update after the one that fills the N-sample window.
        const auto& diag = estimator.diagnostics();
        if(i < N)
        {
            CHECK(diag.used_ekf_fallback);
            continue;
        }
        CHECK_FALSE(diag.used_ekf_fallback);
        CHECK(diag.status == ctrlpp::solve_status::optimal);
    }

    const auto& est = estimator.state();
    CHECK(std::isfinite(est(0)));
    CHECK(std::isfinite(est(1)));
}

namespace
{

constexpr int horizon = 10;
constexpr double eps_mpc = 1e-9;
constexpr std::int32_t tight_budget = 40000;

struct dense_qp
{
    Eigen::MatrixXd P;
    Eigen::VectorXd q;
    Eigen::MatrixXd A;
    Eigen::VectorXd l;
    Eigen::VectorXd u;
};

// The optimum of the posed QP on a working set: the equality rows plus the inequality rows held at a bound. The dual is
// over every constraint row in the backends' sign, P z + q + A^T y = 0, so y > 0 holds an upper bound and y < 0 a lower one.
struct kkt_reference
{
    std::vector<std::pair<Eigen::Index, double>> held;
    Eigen::VectorXd primal;
    Eigen::VectorXd dual;
    Eigen::Index kkt_rows;
    double sigma_min;
    double rounding_bound;
};

struct controller_answer
{
    Eigen::VectorXd primal;
    double first_input;
    ctrlpp::solve_result_status status;
};

auto dense_weight_config(std::optional<double> input_bound) -> ctrlpp::mpc_config<double, NX, NU>
{
    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = horizon,
        .Q = (Eigen::Matrix2d() << 2.0, 0.5, 0.5, 1.0).finished(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
    };
    if(input_bound)
    {
        cfg.u_min = ctrlpp::Vector<double, NU>::Constant(-*input_bound);
        cfg.u_max = ctrlpp::Vector<double, NU>::Constant(*input_bound);
    }
    return cfg;
}

// The QP the regulating solve poses, built by the controller's own builders: Q_state = C^T Q C, the DARE terminal
// weight, and the upper-triangular Hessian expanded to the full symmetric matrix qp_problem's convention denotes.
auto posed_qp(const ctrlpp::discrete_state_space<double, NX, NU, NX>& sys, const ctrlpp::mpc_config<double, NX, NU>& cfg, const Eigen::Vector2d& x0) -> dense_qp
{
    const Eigen::Matrix2d q_state = sys.C.transpose() * cfg.Q * sys.C;
    const auto riccati = ctrlpp::dare<double, NX, NU>(sys.A, sys.B, q_state, cfg.R);
    REQUIRE(riccati.has_value());

    const auto upper = ctrlpp::detail::build_cost_matrix<double, NX, NU>(cfg.horizon, q_state, cfg.R, riccati->P, false, cfg.soft_penalty);
    const Eigen::SparseMatrix<double> full = upper.selfadjointView<Eigen::Upper>();
    const auto A = ctrlpp::detail::build_constraint_matrix<double, NX, NU>(cfg.horizon, sys.A, sys.B, cfg.x_min, cfg.x_max, false, cfg.u_min, cfg.u_max, cfg.du_max);
    auto [l, u] = ctrlpp::detail::build_bounds_vectors<double, NX, NU>(cfg.horizon, x0, cfg.x_min, cfg.x_max, cfg.u_min, cfg.u_max, cfg.du_max);
    return {Eigen::MatrixXd(full), Eigen::VectorXd::Zero(full.cols()), Eigen::MatrixXd(A), std::move(l), std::move(u)};
}

// @cite higham2002 -- Ch. 3 (a computed residual is within gamma_{k+1} (|K| |z| + |rhs|) of the exact one, k the order
// of K) and Ch. 7 (the forward error is at most ||K^{-1}||_2 = 1 / sigma_min(K) times the exact residual).
auto solve_working_set(const dense_qp& qp, std::vector<std::pair<Eigen::Index, double>> held) -> kkt_reference
{
    std::vector<std::pair<Eigen::Index, double>> rows;
    for(Eigen::Index i = 0; i < qp.A.rows(); ++i)
        if(qp.l(i) == qp.u(i))
            rows.emplace_back(i, qp.l(i));
    rows.insert(rows.end(), held.begin(), held.end());

    const Eigen::Index n = qp.P.rows();
    const auto m = static_cast<Eigen::Index>(rows.size());
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(n + m, n + m);
    Eigen::VectorXd rhs(n + m);
    K.topLeftCorner(n, n) = qp.P;
    rhs.head(n) = -qp.q;
    for(Eigen::Index r = 0; r < m; ++r)
    {
        K.block(n + r, 0, 1, n) = qp.A.row(rows[static_cast<std::size_t>(r)].first);
        K.block(0, n + r, n, 1) = qp.A.row(rows[static_cast<std::size_t>(r)].first).transpose();
        rhs(n + r) = rows[static_cast<std::size_t>(r)].second;
    }

    const Eigen::VectorXd z = K.fullPivLu().solve(rhs);
    const double sigma_min = Eigen::JacobiSVD<Eigen::MatrixXd>(K).singularValues().minCoeff();
    const double unit_roundoff = std::numeric_limits<double>::epsilon() / 2.0;
    const double k_terms = static_cast<double>(n + m + 1);
    const double gamma = k_terms * unit_roundoff / (1.0 - k_terms * unit_roundoff);
    const double residual = (K * z - rhs).norm() + gamma * (K.cwiseAbs() * z.cwiseAbs() + rhs.cwiseAbs()).norm();

    Eigen::VectorXd dual = Eigen::VectorXd::Zero(qp.A.rows());
    for(Eigen::Index r = 0; r < m; ++r)
        dual(rows[static_cast<std::size_t>(r)].first) = z(n + r);
    return {std::move(held), z.head(n), std::move(dual), n + m, sigma_min, residual / sigma_min};
}

auto is_inequality(const dense_qp& qp, Eigen::Index row) -> bool
{
    return qp.l(row) < qp.u(row);
}

auto wrong_sign(const kkt_reference& ref, const dense_qp& qp) -> std::optional<std::size_t>
{
    for(std::size_t k = 0; k < ref.held.size(); ++k)
    {
        const auto [row, bound] = ref.held[k];
        const double y = ref.dual(row);
        if((bound == qp.u(row) && y <= 0.0) || (bound == qp.l(row) && y >= 0.0))
            return k;
    }
    return std::nullopt;
}

auto most_violated(const kkt_reference& ref, const dense_qp& qp) -> std::optional<std::pair<Eigen::Index, double>>
{
    std::optional<std::pair<Eigen::Index, double>> worst;
    double worst_gap = 0.0;
    const Eigen::VectorXd Az = qp.A * ref.primal;
    for(Eigen::Index i = 0; i < qp.A.rows(); ++i)
    {
        const bool held = std::ranges::any_of(ref.held, [i](const auto& h) { return h.first == i; });
        if(!is_inequality(qp, i) || held)
            continue;
        const double gap = std::max(qp.l(i) - Az(i), Az(i) - qp.u(i));
        if(gap > worst_gap)
        {
            worst_gap = gap;
            worst = std::pair{i, Az(i) < qp.l(i) ? qp.l(i) : qp.u(i)};
        }
    }
    return worst;
}

// Primal active-set iteration from the unconstrained optimum: release a held bound whose multiplier has the wrong sign,
// else hold the most violated bound. It stops on a working set whose KKT conditions certify the optimum.
auto certified_optimum(const dense_qp& qp) -> kkt_reference
{
    kkt_reference ref = solve_working_set(qp, {});
    const Eigen::Index cycling_guard = 2 * qp.A.rows();
    for(Eigen::Index pass = 0; pass < cycling_guard; ++pass)
    {
        auto held = ref.held;
        if(const auto k = wrong_sign(ref, qp))
            held.erase(held.begin() + static_cast<std::ptrdiff_t>(*k));
        else if(const auto row = most_violated(ref, qp))
            held.push_back(*row);
        else
            return ref;
        ref = solve_working_set(qp, std::move(held));
    }
    FAIL("the active-set iteration did not settle");
    return ref;
}

auto induced_inf_norm(const Eigen::MatrixXd& M) -> double
{
    return M.cwiseAbs().rowwise().sum().maxCoeff();
}

// @cite stellato2020 -- Sec. 3.4: both backends stop once ||A x - z||_inf <= eps_abs + eps_rel max(||A x||, ||z||) and
// ||P x + q + A^T y||_inf <= eps_abs + eps_rel max(||P x||, ||A^T y||, ||q||), on unscaled residuals. On the reference's
// working set those residuals are K delta for the iterate's primal-dual error delta, so ||delta||_2 <= sqrt(m) rho /
// sigma_min with m the KKT rows. Each scale is at most the reference's s plus c ||delta|| (+ rho for z), c the largest
// induced inf-norm of P, A and A^T, so rho <= (eps (1 + s) + eps c ||delta||) / (1 - eps); solving for ||delta||
// gives the bound below. The backends' own floating-point residual evaluation, O(u) relative, is not modeled.
auto stopping_rule_bound(const dense_qp& qp, const kkt_reference& ref, double eps) -> double
{
    const double scale = std::max({(qp.A * ref.primal).lpNorm<Eigen::Infinity>(),
                                   (qp.P * ref.primal).lpNorm<Eigen::Infinity>(),
                                   (qp.A.transpose() * ref.dual).lpNorm<Eigen::Infinity>(),
                                   qp.q.lpNorm<Eigen::Infinity>()});
    const double coupling = std::max({induced_inf_norm(qp.P), induced_inf_norm(qp.A), induced_inf_norm(qp.A.transpose())});
    const double root_m = std::sqrt(static_cast<double>(ref.kkt_rows));
    const double denominator = (1.0 - eps) * ref.sigma_min - root_m * eps * coupling;
    REQUIRE(denominator > 0.0);
    return root_m * eps * (1.0 + scale) / denominator;
}

template <typename Solver>
auto controller_answer_on(const ctrlpp::discrete_state_space<double, NX, NU, NX>& sys, const ctrlpp::mpc_config<double, NX, NU>& cfg, Solver solver, const Eigen::Vector2d& x0) -> controller_answer
{
    auto controller = ctrlpp::test::constructed(ctrlpp::mpc<double, NX, NU, Solver>::create(sys, cfg, std::move(solver)));
    const auto out = controller.solve(x0);
    REQUIRE(out.has_value());
    const auto trajectory = controller.trajectory();
    REQUIRE(trajectory.has_value());

    const auto& [states, inputs] = *trajectory;
    Eigen::VectorXd primal(static_cast<Eigen::Index>(states.size() * NX + inputs.size() * NU));
    Eigen::Index at = 0;
    for(const auto& s : states)
    {
        primal.segment<int(NX)>(at) = s;
        at += static_cast<Eigen::Index>(NX);
    }
    for(const auto& v : inputs)
    {
        primal.segment<int(NU)>(at) = v;
        at += static_cast<Eigen::Index>(NU);
    }
    return {std::move(primal), out->input(0), out->status};
}

void check_against_reference(const controller_answer& answer, const kkt_reference& ref, double bound)
{
    const auto first_input_index = static_cast<Eigen::Index>((horizon + 1) * static_cast<int>(NX));
    INFO("first input " << answer.first_input << ", reference " << ref.primal(first_input_index));
    INFO("primal error " << (answer.primal - ref.primal).norm() << ", bound " << bound + ref.rounding_bound);
    CHECK(answer.status == ctrlpp::solve_result_status::converged);
    CHECK(std::abs(answer.first_input - ref.primal(first_input_index)) <= bound + ref.rounding_bound);
    CHECK((answer.primal - ref.primal).norm() <= bound + ref.rounding_bound);
}

void check_backends_agree(const ctrlpp::mpc_config<double, NX, NU>& cfg, const kkt_reference& ref, const dense_qp& qp, const Eigen::Vector2d& x0)
{
    const auto sys = make_double_integrator();
    const double bound = stopping_rule_bound(qp, ref, eps_mpc);
    INFO("sigma_min " << ref.sigma_min << ", stopping-rule bound " << bound << ", reference rounding bound " << ref.rounding_bound);
    // Polishing is off so each returned primal is the operator-splitting iterate the stopping rule certifies.
    const auto argmin = controller_answer_on(sys, cfg, ctrlpp::argmin_qp_solver{eps_mpc, eps_mpc, tight_budget, false, true, false}, x0);
    check_against_reference(argmin, ref, bound);
#if defined(CTRLPP_HAS_OSQP)
    const auto osqp = controller_answer_on(sys, cfg, ctrlpp::osqp_solver{eps_mpc, eps_mpc, tight_budget, false, true, false}, x0);
    check_against_reference(osqp, ref, bound);
    INFO("argmin to osqp primal distance " << (argmin.primal - osqp.primal).norm());
    CHECK((argmin.primal - osqp.primal).norm() <= 2.0 * bound);
#endif
}

}

TEST_CASE("mpc on argmin_qp_solver agrees with the dense KKT optimum without input bounds", "[mpc][argmin][qp][agree]")
{
    const Eigen::Vector2d x0{1.0, 0.0};
    const auto cfg = dense_weight_config(std::nullopt);
    const auto qp = posed_qp(make_double_integrator(), cfg, x0);
    const auto ref = certified_optimum(qp);
    REQUIRE(ref.held.empty());
    check_backends_agree(cfg, ref, qp, x0);
}

TEST_CASE("mpc on argmin_qp_solver agrees with the KKT-certified optimum with input bounds active", "[mpc][argmin][qp][agree]")
{
    const Eigen::Vector2d x0{1.0, 0.0};
    const auto cfg = dense_weight_config(1.0);
    const auto qp = posed_qp(make_double_integrator(), cfg, x0);
    const auto ref = certified_optimum(qp);
    REQUIRE_FALSE(ref.held.empty());

    // The bound assumes the iterate shares the reference's working set; strict complementarity with every multiplier
    // and slack margin above the bound is the condition for that.
    const double bound = stopping_rule_bound(qp, ref, eps_mpc);
    const Eigen::VectorXd Az = qp.A * ref.primal;
    for(Eigen::Index i = 0; i < qp.A.rows(); ++i)
    {
        if(!is_inequality(qp, i))
            continue;
        const bool held = std::ranges::any_of(ref.held, [i](const auto& h) { return h.first == i; });
        const double margin = held ? std::abs(ref.dual(i)) : std::min(Az(i) - qp.l(i), qp.u(i) - Az(i));
        INFO("row " << i << (held ? " held, multiplier " : " free, slack ") << (held ? ref.dual(i) : margin));
        CHECK(margin > bound + ref.rounding_bound);
    }
    check_backends_agree(cfg, ref, qp, x0);
}

namespace
{

constexpr double eps_qp = 1e-10;

// mhe default-constructs its solver, so the tolerance is bound to the backend type. Polishing is requested, and with no
// box it has no active set to act on.
template <typename Backend, typename SetupError>
class tight_qp_solver
{
public:
    using scalar_type = double;

    tight_qp_solver()
        : inner_{eps_qp, eps_qp, tight_budget, false, true, true}
    {
    }

    auto setup(const ctrlpp::qp_problem<double>& problem) -> ctrlpp::expected<void, SetupError> { return inner_.setup(problem); }

    auto solve(const ctrlpp::qp_update<double>& update) -> ctrlpp::qp_result<double> { return inner_.solve(update); }

private:
    Backend inner_;
};

auto non_diagonal_mhe_config() -> ctrlpp::mhe_config<double, NX, NU, NY, N>
{
    ctrlpp::mhe_config<double, NX, NU, NY, N> cfg;
    cfg.Q = (ctrlpp::Matrix<double, NX, NX>() << 0.02, 0.005, 0.005, 0.01).finished();
    cfg.R = ctrlpp::Matrix<double, NY, NY>::Identity() * 0.1;
    cfg.P0 = ctrlpp::Matrix<double, NX, NX>::Identity() * 10.0;
    return cfg;
}

auto noisy_positions() -> std::vector<ctrlpp::Vector<double, NY>>
{
    std::mt19937 gen(42);
    std::normal_distribution<double> noise(0.0, 0.1);
    ctrlpp::Vector<double, NX> x_true{0.0, 1.0};
    std::vector<ctrlpp::Vector<double, NY>> z;
    for(std::size_t i = 0; i < 3 * N; ++i)
    {
        x_true = linear_dynamics{}(x_true, ctrlpp::Vector<double, NU>::Zero());
        z.push_back((ctrlpp::Vector<double, NY>() << x_true(0) + noise(gen)).finished());
    }
    return z;
}

auto kalman_estimate(const std::vector<ctrlpp::Vector<double, NY>>& z) -> ctrlpp::Vector<double, NX>
{
    const auto cfg = non_diagonal_mhe_config();
    ctrlpp::discrete_state_space<double, NX, NU, NY> sys;
    sys.A << 1.0, dt, 0.0, 1.0;
    sys.B << 0.5 * dt * dt, dt;
    sys.C << 1.0, 0.0;
    sys.D = ctrlpp::Matrix<double, NY, NU>::Zero();
    ctrlpp::kalman_config<double, NX, NU, NY> kf_cfg;
    kf_cfg.Q = cfg.Q;
    kf_cfg.R = cfg.R;
    kf_cfg.P0 = cfg.P0;

    auto filter = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, NX, NU, NY>::create(sys, kf_cfg));
    for(const auto& zk : z)
    {
        filter.predict(ctrlpp::Vector<double, NU>::Zero());
        REQUIRE(filter.update(zk).has_value());
    }
    return filter.state();
}

template <typename Solver>
auto window_estimate(const std::vector<ctrlpp::Vector<double, NY>>& z) -> ctrlpp::Vector<double, NX>
{
    using estimator_t = ctrlpp::mhe<double, NX, NU, NY, N, Solver, linear_dynamics, position_measurement>;
    auto estimator = ctrlpp::test::constructed(estimator_t::create(linear_dynamics{}, position_measurement{}, non_diagonal_mhe_config()));
    for(const auto& zk : z)
    {
        estimator.predict(ctrlpp::Vector<double, NU>::Zero());
        REQUIRE(estimator.update(zk).has_value());
    }
    CHECK_FALSE(estimator.diagnostics().used_ekf_fallback);
    CHECK(estimator.diagnostics().status == ctrlpp::solve_status::optimal);
    return estimator.state();
}

}

TEST_CASE("mhe on argmin_qp_solver agrees with the Kalman filter at the window end", "[mhe][argmin][qp][agree]")
{
    const auto z = noisy_positions();
    const auto reference = kalman_estimate(z);

    // The window-end marginal of the unconstrained window QP is the Kalman posterior; a splitting solver stopped at
    // eps_qp reaches it to eps_qp times the solution scale, propagated through the N coupled stage blocks.
    const double tol = static_cast<double>(N) * eps_qp * (1.0 + reference.norm());
    const auto argmin = window_estimate<tight_qp_solver<ctrlpp::argmin_qp_solver, ctrlpp::argmin_qp_setup_error>>(z);
    INFO("argmin error " << (argmin - reference).cwiseAbs().maxCoeff() << ", bound " << tol);
    CHECK((argmin - reference).cwiseAbs().maxCoeff() <= tol);
#if defined(CTRLPP_HAS_OSQP)
    const auto osqp = window_estimate<tight_qp_solver<ctrlpp::osqp_solver, ctrlpp::osqp_setup_error>>(z);
    INFO("osqp error " << (osqp - reference).cwiseAbs().maxCoeff() << ", argmin to osqp " << (argmin - osqp).cwiseAbs().maxCoeff());
    CHECK((osqp - reference).cwiseAbs().maxCoeff() <= tol);
    CHECK((argmin - osqp).cwiseAbs().maxCoeff() <= 2.0 * tol);
#endif
}

#endif
