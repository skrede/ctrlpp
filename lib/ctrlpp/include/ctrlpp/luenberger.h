#ifndef HPP_GUARD_CTRLPP_LUENBERGER_H
#define HPP_GUARD_CTRLPP_LUENBERGER_H

#include "ctrlpp/types.h"
#include "ctrlpp/state_space.h"
#include "ctrlpp/observer_policy.h"

#include <cstddef>
#include <utility>

namespace ctrlpp {

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
class luenberger_observer {
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
        : sys_{std::move(sys)}
        , L_{std::move(L)}
        , x_{std::move(x0)}
    {
    }

    void predict(const input_vector_t& u)
    {
        x_ = (sys_.A * x_ + sys_.B * u).eval();
    }

    void update(const output_vector_t& z)
    {
        x_ = (x_ + L_ * (z - sys_.C * x_)).eval();
    }

    [[nodiscard]] auto state() const -> const state_vector_t& { return x_; }

    void set_gain(const gain_matrix_t& L) { L_ = L; }
    void set_model(system_t sys) { sys_ = std::move(sys); }
    void reset(const state_vector_t& x0) { x_ = x0; }

private:
    system_t sys_;
    gain_matrix_t L_;
    state_vector_t x_;
};

static_assert(ObserverPolicy<luenberger_observer<double, 2, 1, 1>>);
static_assert(!CovarianceObserver<luenberger_observer<double, 2, 1, 1>>);

}

#endif
