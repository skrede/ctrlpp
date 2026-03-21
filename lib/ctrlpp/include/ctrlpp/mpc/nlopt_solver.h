#ifndef HPP_GUARD_CTRLPP_MPC_NLOPT_SOLVER_H
#define HPP_GUARD_CTRLPP_MPC_NLOPT_SOLVER_H

#include "ctrlpp/mpc/nlp_solver.h"

#include <nlopt.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace ctrlpp {

/// NLopt algorithm selection for the NMPC solver backend.
/// Note: MMA does not support equality constraints.
enum class nlopt_algorithm : std::uint8_t {
    slsqp,
    mma,
    cobyla,
    isres
};

template<typename Scalar>
struct nlopt_settings {
    nlopt_algorithm algorithm{nlopt_algorithm::slsqp};
    Scalar ftol_rel{Scalar{1e-6}};
    Scalar xtol_rel{Scalar{1e-6}};
    int max_eval{500};
    Scalar max_time{Scalar{0}};
    Scalar constraint_tol{Scalar{1e-8}};
};

template<typename Scalar>
class nlopt_solver {
    static_assert(std::is_same_v<Scalar, double>,
                  "NLopt operates in double precision only");

public:
    using scalar_type = Scalar;

    explicit nlopt_solver(nlopt_settings<Scalar> settings = {})
        : settings_{settings}
    {
    }

    void setup(const nlp_problem<Scalar>& problem) {
        problem_ = &problem;

        opt_ = nlopt::opt(to_nlopt_algorithm(settings_.algorithm), problem.n_vars);

        // Variable bounds
        if (problem.x_lower.size() > 0) {
            opt_.set_lower_bounds(to_stdvec(problem.x_lower));
        }
        if (problem.x_upper.size() > 0) {
            opt_.set_upper_bounds(to_stdvec(problem.x_upper));
        }

        // Objective
        opt_.set_min_objective(objective_callback, this);

        // Partition constraints into equality and inequality
        if (problem.n_constraints > 0) {
            eq_indices_.clear();
            ineq_lower_indices_.clear();
            ineq_upper_indices_.clear();

            for (int i = 0; i < problem.n_constraints; ++i) {
                if (problem.c_lower[i] == problem.c_upper[i]) {
                    eq_indices_.push_back(i);
                } else {
                    if (std::isfinite(static_cast<double>(problem.c_upper[i]))) {
                        ineq_upper_indices_.push_back(i);
                    }
                    if (std::isfinite(static_cast<double>(problem.c_lower[i]))) {
                        ineq_lower_indices_.push_back(i);
                    }
                }
            }

            if (settings_.algorithm == nlopt_algorithm::mma && !eq_indices_.empty()) {
                throw std::invalid_argument(
                    "MMA algorithm does not support equality constraints");
            }

            // Equality constraints: c(x) - target = 0
            if (!eq_indices_.empty()) {
                std::vector<double> eq_tol(eq_indices_.size(),
                                           static_cast<double>(settings_.constraint_tol));
                opt_.add_equality_mconstraint(equality_callback, this, eq_tol);
            }

            // Inequality constraints: c(x) - c_upper <= 0
            if (!ineq_upper_indices_.empty()) {
                std::vector<double> tol(ineq_upper_indices_.size(),
                                        static_cast<double>(settings_.constraint_tol));
                opt_.add_inequality_mconstraint(ineq_upper_callback, this, tol);
            }

            // Inequality constraints: c_lower - c(x) <= 0
            if (!ineq_lower_indices_.empty()) {
                std::vector<double> tol(ineq_lower_indices_.size(),
                                        static_cast<double>(settings_.constraint_tol));
                opt_.add_inequality_mconstraint(ineq_lower_callback, this, tol);
            }
        }

        // Stopping criteria
        opt_.set_ftol_rel(static_cast<double>(settings_.ftol_rel));
        opt_.set_xtol_rel(static_cast<double>(settings_.xtol_rel));
        opt_.set_maxeval(settings_.max_eval);
        if (settings_.max_time > Scalar{0}) {
            opt_.set_maxtime(static_cast<double>(settings_.max_time));
        }
    }

    auto solve(const nlp_update<Scalar>& update) -> nlp_result<Scalar> {
        nlp_result<Scalar> result{};
        auto x = to_stdvec(update.x0);
        double opt_f{};

        auto t0 = std::chrono::steady_clock::now();

        try {
            auto rc = opt_.optimize(x, opt_f);
            result.status = map_result(rc);
        } catch (const nlopt::roundoff_limited&) {
            result.status = solve_status::solved_inaccurate;
            opt_f = problem_->cost({x.data(), x.size()});
        } catch (const std::runtime_error&) {
            result.status = solve_status::error;
        }

        auto t1 = std::chrono::steady_clock::now();
        result.solve_time = static_cast<Scalar>(
            std::chrono::duration<double>(t1 - t0).count());

        // Copy solution back to Eigen
        result.x.resize(static_cast<Eigen::Index>(x.size()));
        for (Eigen::Index i = 0; i < result.x.size(); ++i) {
            result.x[i] = static_cast<Scalar>(x[static_cast<std::size_t>(i)]);
        }

        result.objective = static_cast<Scalar>(opt_f);
        result.iterations = opt_.get_numevals();

        // Compute primal residual as max constraint violation
        result.primal_residual = compute_max_violation(x);

        return result;
    }

private:
    static auto to_nlopt_algorithm(nlopt_algorithm alg) -> nlopt::algorithm {
        switch (alg) {
        case nlopt_algorithm::slsqp: return nlopt::LD_SLSQP;
        case nlopt_algorithm::mma:   return nlopt::LD_MMA;
        case nlopt_algorithm::cobyla: return nlopt::LN_COBYLA;
        case nlopt_algorithm::isres: return nlopt::GN_ISRES;
        }
        return nlopt::LD_SLSQP;
    }

    static auto to_stdvec(const Eigen::VectorX<Scalar>& v) -> std::vector<double> {
        std::vector<double> out(static_cast<std::size_t>(v.size()));
        for (Eigen::Index i = 0; i < v.size(); ++i) {
            out[static_cast<std::size_t>(i)] = static_cast<double>(v[i]);
        }
        return out;
    }

    static constexpr auto map_result(nlopt::result rc) -> solve_status {
        switch (rc) {
        case nlopt::SUCCESS:          [[fallthrough]];
        case nlopt::STOPVAL_REACHED:  [[fallthrough]];
        case nlopt::FTOL_REACHED:     [[fallthrough]];
        case nlopt::XTOL_REACHED:     return solve_status::optimal;
        case nlopt::MAXEVAL_REACHED:  return solve_status::max_iterations;
        case nlopt::MAXTIME_REACHED:  return solve_status::time_limit;
        case nlopt::ROUNDOFF_LIMITED: return solve_status::solved_inaccurate;
        default:                      return solve_status::error;
        }
    }

    static auto objective_callback(const std::vector<double>& x,
                                   std::vector<double>& grad,
                                   void* data) -> double {
        auto* self = static_cast<nlopt_solver*>(data);
        std::span<const double> x_span{x.data(), x.size()};

        if (!grad.empty() && self->problem_->gradient) {
            std::span<double> grad_span{grad.data(), grad.size()};
            self->problem_->gradient(x_span, grad_span);
        }

        return self->problem_->cost(x_span);
    }

    static void equality_callback(unsigned m, double* result,
                                  unsigned n, const double* x,
                                  double* grad, void* data) {
        auto* self = static_cast<nlopt_solver*>(data);

        // Evaluate all constraints
        self->eval_constraints(n, x);

        // Extract equality constraints: c(x) - target = 0
        for (unsigned i = 0; i < m; ++i) {
            auto idx = self->eq_indices_[i];
            result[i] = self->constraint_buf_[static_cast<std::size_t>(idx)]
                      - static_cast<double>(self->problem_->c_lower[idx]);
        }

        if (grad) {
            self->fd_mconstraint_grad(m, n, x, grad, self->eq_indices_,
                [](double ci, double bound) { return ci - bound; },
                [self](int idx) {
                    return static_cast<double>(self->problem_->c_lower[idx]);
                });
        }
    }

    static void ineq_upper_callback(unsigned m, double* result,
                                    unsigned n, const double* x,
                                    double* grad, void* data) {
        auto* self = static_cast<nlopt_solver*>(data);
        self->eval_constraints(n, x);

        // c(x) - c_upper <= 0
        for (unsigned i = 0; i < m; ++i) {
            auto idx = self->ineq_upper_indices_[i];
            result[i] = self->constraint_buf_[static_cast<std::size_t>(idx)]
                      - static_cast<double>(self->problem_->c_upper[idx]);
        }

        if (grad) {
            self->fd_mconstraint_grad(m, n, x, grad, self->ineq_upper_indices_,
                [](double ci, double bound) { return ci - bound; },
                [self](int idx) {
                    return static_cast<double>(self->problem_->c_upper[idx]);
                });
        }
    }

    static void ineq_lower_callback(unsigned m, double* result,
                                    unsigned n, const double* x,
                                    double* grad, void* data) {
        auto* self = static_cast<nlopt_solver*>(data);
        self->eval_constraints(n, x);

        // c_lower - c(x) <= 0
        for (unsigned i = 0; i < m; ++i) {
            auto idx = self->ineq_lower_indices_[i];
            result[i] = static_cast<double>(self->problem_->c_lower[idx])
                      - self->constraint_buf_[static_cast<std::size_t>(idx)];
        }

        if (grad) {
            self->fd_mconstraint_grad(m, n, x, grad, self->ineq_lower_indices_,
                [](double ci, double bound) { return bound - ci; },
                [self](int idx) {
                    return static_cast<double>(self->problem_->c_lower[idx]);
                });
        }
    }

    template<typename TransformFn, typename BoundFn>
    void fd_mconstraint_grad(unsigned m, unsigned n, const double* x,
                             double* grad,
                             const std::vector<int>& indices,
                             TransformFn transform,
                             BoundFn get_bound) {
        const auto eps = std::sqrt(std::numeric_limits<double>::epsilon());
        fd_x_buf_.assign(x, x + n);
        fd_c_plus_.resize(static_cast<std::size_t>(problem_->n_constraints));
        fd_c_minus_.resize(static_cast<std::size_t>(problem_->n_constraints));

        for (unsigned j = 0; j < n; ++j) {
            const double h = eps * std::max(1.0, std::abs(x[j]));
            const double orig = fd_x_buf_[j];

            fd_x_buf_[j] = orig + h;
            std::span<const double> xp{fd_x_buf_.data(), n};
            std::span<double> cp{fd_c_plus_.data(), fd_c_plus_.size()};
            problem_->constraints(xp, cp);

            fd_x_buf_[j] = orig - h;
            std::span<const double> xm{fd_x_buf_.data(), n};
            std::span<double> cm{fd_c_minus_.data(), fd_c_minus_.size()};
            problem_->constraints(xm, cm);

            for (unsigned i = 0; i < m; ++i) {
                auto idx = indices[i];
                double bound = get_bound(idx);
                double fp = transform(fd_c_plus_[static_cast<std::size_t>(idx)], bound);
                double fm = transform(fd_c_minus_[static_cast<std::size_t>(idx)], bound);
                grad[static_cast<std::size_t>(i) * n + j] = (fp - fm) / (2.0 * h);
            }

            fd_x_buf_[j] = orig;
        }
    }

    void eval_constraints(unsigned n, const double* x) {
        constraint_buf_.resize(static_cast<std::size_t>(problem_->n_constraints));
        std::span<const double> x_span{x, n};
        std::span<double> c_span{constraint_buf_.data(), constraint_buf_.size()};
        problem_->constraints(x_span, c_span);
    }

    auto compute_max_violation(const std::vector<double>& x) const -> Scalar {
        if (problem_->n_constraints == 0) {
            return Scalar{0};
        }

        eval_constraints_const(x);

        Scalar max_viol{0};
        for (int i = 0; i < problem_->n_constraints; ++i) {
            auto ci = static_cast<Scalar>(
                constraint_buf_[static_cast<std::size_t>(i)]);
            Scalar viol{0};
            if (ci > problem_->c_upper[i]) {
                viol = ci - problem_->c_upper[i];
            } else if (ci < problem_->c_lower[i]) {
                viol = problem_->c_lower[i] - ci;
            }
            max_viol = std::max(max_viol, viol);
        }
        return max_viol;
    }

    void eval_constraints_const(const std::vector<double>& x) const {
        constraint_buf_.resize(static_cast<std::size_t>(problem_->n_constraints));
        std::span<const double> x_span{x.data(), x.size()};
        std::span<double> c_span{constraint_buf_.data(), constraint_buf_.size()};
        problem_->constraints(x_span, c_span);
    }

    nlopt_settings<Scalar> settings_;
    const nlp_problem<Scalar>* problem_{nullptr};
    nlopt::opt opt_{};

    std::vector<int> eq_indices_;
    std::vector<int> ineq_upper_indices_;
    std::vector<int> ineq_lower_indices_;
    mutable std::vector<double> constraint_buf_;
    mutable std::vector<double> fd_x_buf_;
    mutable std::vector<double> fd_c_plus_;
    mutable std::vector<double> fd_c_minus_;
};

}

#endif
