#ifndef HPP_GUARD_CTRLPP_SYSID_RECURSIVE_ARX_H
#define HPP_GUARD_CTRLPP_SYSID_RECURSIVE_ARX_H

#include "ctrlpp/sysid/rls.h"
#include "ctrlpp/state_space.h"
#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <array>
#include <cstddef>

namespace ctrlpp {

template<typename Scalar, std::size_t NA, std::size_t NB, std::size_t NU = 1, std::size_t NY = 1>
class recursive_arx {
  public:
    static constexpr std::size_t NP = NA * NY + NB * NU;

    explicit recursive_arx(rls_config<Scalar, NP> config = {})
        : rls_{config}
    {
    }

    void update(Scalar y, Scalar u)
    {
        Vector<Scalar, NP> phi = Vector<Scalar, NP>::Zero();

        // Build regressor: [y(t-1), ..., y(t-NA), u(t-1), ..., u(t-NB)]
        for (std::size_t i = 0; i < NA; ++i) {
            std::size_t idx = (write_idx_ + NA - 1 - i) % NA;
            phi(static_cast<int>(i)) = y_hist_[idx];
        }
        for (std::size_t i = 0; i < NB; ++i) {
            std::size_t idx = (write_idx_ + NB - 1 - i) % NB;
            phi(static_cast<int>(NA + i)) = u_hist_[idx];
        }

        rls_.update(y, phi);

        y_hist_[write_idx_ % NA] = y;
        u_hist_[write_idx_ % NB] = u;
        ++write_idx_;
        ++sample_count_;
    }

    [[nodiscard]] auto parameters() const -> const Vector<Scalar, NP>&
    {
        return rls_.parameters();
    }

    [[nodiscard]] auto covariance() const -> const Matrix<Scalar, NP, NP>&
    {
        return rls_.covariance();
    }

    [[nodiscard]] auto to_state_space() const -> discrete_state_space<Scalar, NA, NU, NY>
    {
        auto theta = rls_.parameters();

        Matrix<Scalar, NA, NA> A = Matrix<Scalar, NA, NA>::Zero();
        Matrix<Scalar, NA, NU> B = Matrix<Scalar, NA, NU>::Zero();
        Matrix<Scalar, NY, NA> C = Matrix<Scalar, NY, NA>::Zero();
        Matrix<Scalar, NY, NU> D = Matrix<Scalar, NY, NU>::Zero();

        // Observer canonical form for ARX:
        //   y(t) = a1*y(t-1) + ... + aNa*y(t-NA) + b1*u(t-1) + ... + bNb*u(t-NB)
        //
        //   A = [a1  1  0 ...]    B = [b1]    C = [1 0 ... 0]    D = [0]
        //       [a2  0  1 ...]        [b2]
        //       [... ... ...]         [...]
        //       [aNa 0  0 ...]        [bNa (or 0)]
        //
        // A: first column = a-coefficients, superdiagonal = 1
        for (std::size_t i = 0; i < NA; ++i) {
            A(static_cast<int>(i), 0) = theta(static_cast<int>(i));
        }
        for (std::size_t i = 0; i + 1 < NA; ++i) {
            A(static_cast<int>(i), static_cast<int>(i + 1)) = Scalar{1};
        }

        // B column: b-coefficients, zero-padded if NB < NA
        for (std::size_t i = 0; i < NB && i < NA; ++i) {
            B(static_cast<int>(i), 0) = theta(static_cast<int>(NA + i));
        }

        // C = [1, 0, ..., 0]
        C(0, 0) = Scalar{1};

        return {.A = A, .B = B, .C = C, .D = D};
    }

  private:
    rls<Scalar, NP> rls_;
    std::array<Scalar, NA> y_hist_{};
    std::array<Scalar, NB> u_hist_{};
    std::size_t write_idx_{0};
    std::size_t sample_count_{0};
};

}

#endif
