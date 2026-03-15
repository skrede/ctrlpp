#ifndef HPP_GUARD_CPPCTRL_SOLVER_POLICY_H
#define HPP_GUARD_CPPCTRL_SOLVER_POLICY_H

namespace ctrlpp {

template<typename P>
concept SolverPolicy = requires {
    typename P::solver_tag;
};

}

#endif
