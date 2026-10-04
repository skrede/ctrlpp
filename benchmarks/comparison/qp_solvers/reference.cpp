#include "judge.h"
#include "reference.h"
#include "argmin_arms.h"
#include "competitor_arms.h"

#include <cmath>
#include <memory>
#include <algorithm>
#include <functional>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

namespace qp = problems::qp;

/// One reference solver on one form, posed on first use and warm-started from its own previous tick.
struct candidate
{
    std::unique_ptr<qp_arm> arm;
    std::reference_wrapper<const problem_case> form;
    bool posed;
    qp::dense_program program;
    qp_tick tick;
};

auto make_candidate(std::unique_ptr<qp_arm> arm, const problem_case& form) -> candidate
{
    const qp_instance& instance = form.instance();
    return {std::move(arm), form, false, {instance.P, {}, instance.A, {}, {}}, {}};
}

auto try_candidate(candidate& c, const Eigen::VectorXd& x0, reference_step& out) -> bool
{
    const problem_case& form = c.form.get();
    form.fill(x0, c.tick);
    if(!c.posed)
        c.arm->setup(form.instance(), c.tick);
    c.posed = true;
    c.arm->update(c.tick);
    c.arm->solve();
    arm_answer answer;
    c.arm->answer(answer);
    c.program.q = c.tick.q;
    c.program.l = c.tick.l;
    c.program.u = c.tick.u;
    if(!(judge(c.program, answer.x, answer.y).worst <= reference_kkt_max))
        return false;
    out = {true, row_classes(form.rows_of(c.tick, answer.x))};
    return true;
}

auto candidates_for(const reference_forms& forms, hpipm_dual_signs signs) -> std::vector<candidate>
{
    std::vector<candidate> list;
    list.push_back(make_candidate(make_qpoases_arm(tier::high), forms.dense));
    if(forms.sparse != nullptr)
        list.push_back(make_candidate(make_hpipm_reference(signs), *forms.sparse));
    list.push_back(make_candidate(make_argmin_active_set_reference(), forms.dense));
    return list;
}

}

auto row_classes(const inequality_rows& rows) -> std::vector<qp::row_state>
{
    std::vector<qp::row_state> classes(static_cast<std::size_t>(rows.value.size()), qp::row_state::inactive);
    for(Eigen::Index i = 0; i < rows.value.size(); ++i)
    {
        const double scale = 1.0 + std::max(std::abs(rows.lower(i)), std::abs(rows.upper(i)));
        const double to_lower = rows.value(i) - rows.lower(i);
        const double to_upper = rows.upper(i) - rows.value(i);
        if(std::min(to_lower, to_upper) <= binding_row_tolerance * scale)
            classes[static_cast<std::size_t>(i)] = to_lower <= to_upper ? qp::row_state::at_lower : qp::row_state::at_upper;
    }
    return classes;
}

auto reference_sequence(const reference_forms& forms, regime kind, std::int32_t seed, std::int32_t steps,
                        hpipm_dual_signs signs) -> std::vector<reference_step>
{
    std::vector<candidate> list = candidates_for(forms, signs);
    std::vector<reference_step> sequence(static_cast<std::size_t>(steps), reference_step{false, {}});
    const std::int32_t distinct = kind == regime::fixed ? std::min<std::int32_t>(steps, 1) : steps;
    for(std::int32_t s = 0; s < distinct; ++s)
    {
        const Eigen::VectorXd x0 = initial_state(forms.dense, kind, seed, s);
        for(candidate& c : list)
            if(try_candidate(c, x0, sequence[static_cast<std::size_t>(s)]))
                break;
    }
    std::fill(sequence.begin() + distinct, sequence.end(), distinct > 0 ? sequence.front() : reference_step{false, {}});
    return sequence;
}

auto binding_fraction(const reference_step& step) -> double
{
    if(!step.found)
        return -1.0;
    if(step.rows.empty())
        return 0.0;
    const auto binding = std::ranges::count_if(step.rows, [](qp::row_state s) { return s != qp::row_state::inactive; });
    return static_cast<double>(binding) / static_cast<double>(step.rows.size());
}

auto churn_fraction(const reference_step& previous, const reference_step& current) -> double
{
    if(!previous.found || !current.found || current.rows.empty())
        return -1.0;
    std::size_t changed = 0;
    for(std::size_t i = 0; i < current.rows.size(); ++i)
        changed += previous.rows[i] != current.rows[i] ? 1 : 0;
    return static_cast<double>(changed) / static_cast<double>(current.rows.size());
}

}
