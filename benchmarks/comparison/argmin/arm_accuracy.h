#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_ARGMIN_ARM_ACCURACY_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_ARGMIN_ARM_ACCURACY_H

#include <Eigen/Core>

#include <limits>
#include <string>
#include <vector>
#include <cstddef>
#include <algorithm>

namespace ctrlpp::bench::argmin_arms
{

constexpr char const* spread_metric = "max abs entrywise spread of the arms' solutions";
constexpr char const* deviation_metric = "max abs entrywise deviation of the two arms' solutions";
constexpr char const* violation_metric = "max abs constraint violation of this arm's own returned solution";

/// One arm's answer, held so a set of arms is scored against each other only
/// after every one of them has produced one.
struct arm_answer
{
    double violation;
    std::string label;
    Eigen::VectorXd solution;
};

/// A pairwise deviation is not expressible over more than two arms, so a set
/// publishes the largest entrywise spread of its solutions. With two arms that
/// is exactly the deviation, which is the test that this generalizes the
/// two-arm figure rather than replacing it.
///
/// An arm that exits on a budget leaves no solution behind, and dropping it
/// from the spread would publish an agreement figure for a set that never
/// agreed. The spread is reported as not-a-number instead.
inline auto solution_spread(const std::vector<arm_answer>& arms) -> double
{
    if(arms.empty())
        return std::numeric_limits<double>::quiet_NaN();
    for(const arm_answer& arm : arms)
        if(arm.solution.size() == 0 || arm.solution.size() != arms.front().solution.size())
            return std::numeric_limits<double>::quiet_NaN();

    double spread = 0.0;
    for(std::size_t i = 0; i < arms.size(); ++i)
        for(std::size_t j = i + 1; j < arms.size(); ++j)
            spread = std::max(spread, (arms[i].solution - arms[j].solution).cwiseAbs().maxCoeff());
    return spread;
}

inline auto answer_for_label(const std::vector<arm_answer>& arms, const std::string& label) -> const arm_answer*
{
    for(const arm_answer& arm : arms)
        if(arm.label == label)
            return &arm;
    return nullptr;
}

}

#endif
