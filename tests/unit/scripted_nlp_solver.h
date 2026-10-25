#ifndef HPP_GUARD_CTRLPP_TESTS_SCRIPTED_NLP_SOLVER_H
#define HPP_GUARD_CTRLPP_TESTS_SCRIPTED_NLP_SOLVER_H

// Backend-free NLP solver whose answer is scripted call by call, for testing what
// a nonlinear controller commits after each answer. The script and everything the
// backend saw live in a shared record, so they stay observable after the solver
// has been moved into a controller.
//
// Each call takes the next script entry; the last entry repeats once the script is
// exhausted, and an empty script answers optimal with nothing corrupted. The
// finite primal of call c at index i is (c + 1) * (i + 1) / 64, exact in binary
// and distinct per call and per index, so a warm-start shift or an extracted
// input shows which call's primal it came from. With SolveInto set, the solver
// also offers the write-into shape, which leaves nlp_result::x empty.

#include "ctrlpp/mpc/nlp_solver.h"

#include "stub_qp_solver.h"

#include <Eigen/Core>

#include <limits>
#include <memory>
#include <vector>
#include <cstddef>
#include <utility>
#include <algorithm>

namespace ctrlpp_test
{

enum class corrupted_field
{
    none,
    primal,
    objective,
    solve_time,
    primal_residual,
    iterations,
    all
};

struct scripted_answer
{
    // A negative or positive length offset makes the by-value primal shorter or
    // longer than the posed problem; the write-into shape cannot change length.
    explicit scripted_answer(ctrlpp::solve_status reported_status,
                             corrupted_field corrupted = corrupted_field::none,
                             report_values corrupting = report_values::nan,
                             int offset = 0)
        : status{reported_status}
        , field{corrupted}
        , value{corrupting}
        , length_offset{offset}
    {
    }

    ctrlpp::solve_status status;
    corrupted_field field;
    report_values value;
    int length_offset;
};

template <typename Scalar>
struct scripted_record
{
    std::vector<scripted_answer> script;
    std::vector<Eigen::VectorX<Scalar>> warm_starts;
    std::vector<Eigen::VectorX<Scalar>> primals;
};

template <typename Scalar>
auto corrupting_value(report_values value) -> Scalar
{
    switch(value)
    {
    case report_values::positive_infinity:
        return std::numeric_limits<Scalar>::infinity();
    case report_values::negative_infinity:
        return -std::numeric_limits<Scalar>::infinity();
    case report_values::finite:
        return Scalar{1};
    case report_values::nan:
    default:
        return std::numeric_limits<Scalar>::quiet_NaN();
    }
}

template <typename Scalar, bool SolveInto = false>
class scripted_nlp_solver
{
public:
    using scalar_type = Scalar;

    scripted_nlp_solver()
        : m_record{std::make_shared<scripted_record<Scalar>>()}
        , m_n_vars{0}
    {
    }

    explicit scripted_nlp_solver(std::shared_ptr<scripted_record<Scalar>> record)
        : m_record{std::move(record)}
        , m_n_vars{0}
    {
    }

    auto setup(const ctrlpp::nlp_problem<Scalar>& problem) -> ctrlpp::expected<void, stub_setup_error>
    {
        m_n_vars = problem.n_vars;
        return {};
    }

    template <int NV>
    auto setup(const ctrlpp::nlp_problem_static<Scalar, NV>& problem) -> ctrlpp::expected<void, stub_setup_error>
    {
        m_n_vars = problem.n_vars;
        return {};
    }

    auto solve(const ctrlpp::nlp_update<Scalar>& update) -> ctrlpp::nlp_result<Scalar>
    {
        const scripted_answer answer = next_answer(update);
        ctrlpp::nlp_result<Scalar> result = report(answer);
        result.x = primal(answer, std::max(m_n_vars + answer.length_offset, 0));
        m_record->primals.push_back(result.x);
        return result;
    }

    auto solve_into(const ctrlpp::nlp_update<Scalar>& update, Eigen::Ref<Eigen::VectorX<Scalar>> solution) -> ctrlpp::nlp_result<Scalar>
        requires SolveInto
    {
        const scripted_answer answer = next_answer(update);
        solution = primal(answer, static_cast<int>(solution.size()));
        m_record->primals.push_back(solution);
        return report(answer);
    }

private:
    std::shared_ptr<scripted_record<Scalar>> m_record;
    int m_n_vars;

    auto next_answer(const ctrlpp::nlp_update<Scalar>& update) -> scripted_answer
    {
        const std::size_t call = m_record->warm_starts.size();
        m_record->warm_starts.push_back(update.x0);
        const auto& script = m_record->script;
        if(script.empty())
            return scripted_answer{ctrlpp::solve_status::optimal};
        return script[std::min(call, script.size() - 1)];
    }

    auto primal(const scripted_answer& answer, int length) const -> Eigen::VectorX<Scalar>
    {
        const auto call = static_cast<Scalar>(m_record->warm_starts.size());
        Eigen::VectorX<Scalar> x = Eigen::VectorX<Scalar>::LinSpaced(length, Scalar{1}, static_cast<Scalar>(length)) * (call / Scalar{64});
        if(corrupts(answer, corrupted_field::primal))
            x.setConstant(corrupting_value<Scalar>(answer.value));
        return x;
    }

    static auto corrupts(const scripted_answer& answer, corrupted_field field) -> bool { return answer.field == field || answer.field == corrupted_field::all; }

    static auto report(const scripted_answer& answer) -> ctrlpp::nlp_result<Scalar>
    {
        const Scalar bad = corrupting_value<Scalar>(answer.value);
        ctrlpp::nlp_result<Scalar> result{};
        result.status = answer.status;
        result.objective = corrupts(answer, corrupted_field::objective) ? bad : Scalar{1};
        result.solve_time = corrupts(answer, corrupted_field::solve_time) ? bad : Scalar{0};
        result.iterations = corrupts(answer, corrupted_field::iterations) ? -1 : 1;
        result.primal_residual = corrupts(answer, corrupted_field::primal_residual) ? bad : Scalar{0};
        return result;
    }
};

}

#endif
