#ifndef HPP_GUARD_CPPCTRL_OBSERVER_POLICY_H
#define HPP_GUARD_CPPCTRL_OBSERVER_POLICY_H

#include <concepts>
#include <variant>

namespace ctrlpp {

template<typename O>
concept ObserverPolicy = requires {
    typename O::observer_tag;
    typename O::state_vector_t;
    typename O::input_vector_t;
    typename O::output_vector_t;
} && requires(O obs,
              const typename O::input_vector_t& u,
              const typename O::output_vector_t& y) {
    obs.predict(u);
    obs.update(y);
    { obs.state() } -> std::convertible_to<const typename O::state_vector_t&>;
};

template<typename O>
concept CovarianceObserver = ObserverPolicy<O> && requires(const O& obs) {
    { obs.covariance() };
    { obs.innovation() };
};

struct NullObserver {
    using observer_tag = void;
    using state_vector_t = std::monostate;
    using input_vector_t = std::monostate;
    using output_vector_t = std::monostate;

    void predict(const std::monostate&) {}
    void update(const std::monostate&) {}
    auto state() const -> const std::monostate& { static std::monostate s; return s; }
};

}

#endif
