#ifndef HPP_GUARD_CTRLPP_MODEL_DISCRETIZE_H
#define HPP_GUARD_CTRLPP_MODEL_DISCRETIZE_H

// @cite franklin2015 -- Franklin et al., "Feedback Control of Dynamic Systems", 2015
// @cite astrom1997 -- Astrom & Wittenmark, "Computer-Controlled Systems: Theory and Design", 3rd ed., 1997, Ch. 3

#include "ctrlpp/model/state_space.h"

#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <unsupported/Eigen/MatrixFunctions>

#include <cmath>
#include <cstddef>

namespace ctrlpp
{

struct zoh
{
};

struct tustin
{
};

template <typename Scalar>
struct tustin_prewarp
{
    Scalar w_c; // rad/s
};

struct forward_euler
{
};

struct backward_euler
{
};

// Van Loan's method: Ad and Bd are the upper blocks of exp([[A*dt, B*dt], [0, 0]]).
// @cite vanloan1978 -- Van Loan, "Computing Integrals Involving the Matrix Exponential",
//   IEEE Trans. Autom. Control 23(3):395-404, 1978
// @cite astrom1997 -- Astrom & Wittenmark, "Computer-Controlled Systems", 3rd ed., 1997, Sec. 3.2
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY> discretize(zoh, const continuous_state_space<Scalar, NX, NU, NY>& sys, Scalar dt)
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);
    constexpr int aug = nx + nu;

    Eigen::Matrix<Scalar, aug, aug> M;
    M.setZero();
    M.template block<nx, nx>(0, 0) = sys.A * dt;
    M.template block<nx, nu>(0, nx) = sys.B * dt;

    Eigen::Matrix<Scalar, aug, aug> expM = M.exp();

    Matrix<Scalar, NX, NX> Ad = expM.template block<nx, nx>(0, 0);
    Matrix<Scalar, NX, NU> Bd = expM.template block<nx, nu>(0, nx);

    return {Ad, Bd, sys.C, sys.D};
}

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY> discretize(const continuous_state_space<Scalar, NX, NU, NY>& sys, Scalar dt, zoh = {})
{
    return discretize(zoh{}, sys, dt);
}

namespace detail
{

// Dd = D + C*Bd/2 is the feed-through a bilinear map adds even to a strictly proper continuous system.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY> tustin_bilinear(const continuous_state_space<Scalar, NX, NU, NY>& sys, Scalar dt)
{
    const Matrix<Scalar, NX, NX> identity = Matrix<Scalar, NX, NX>::Identity();
    const Matrix<Scalar, NX, NX> half_a_dt = sys.A * (dt / Scalar{2});
    const Matrix<Scalar, NX, NX> forward_resolvent_inverse = (identity - half_a_dt).inverse();

    Matrix<Scalar, NX, NX> Ad = forward_resolvent_inverse * (identity + half_a_dt);
    Matrix<Scalar, NX, NU> Bd = forward_resolvent_inverse * sys.B * dt;
    Matrix<Scalar, NY, NX> Cd = sys.C * forward_resolvent_inverse;
    Matrix<Scalar, NY, NU> Dd = sys.D + Cd * sys.B * (dt / Scalar{2});

    return {Ad, Bd, Cd, Dd};
}

}

// @cite astrom1997 -- Astrom & Wittenmark, "Computer-Controlled Systems", 3rd ed., 1997, Sec. 3.5
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY> discretize(tustin, const continuous_state_space<Scalar, NX, NU, NY>& sys, Scalar dt)
{
    return detail::tustin_bilinear(sys, dt);
}

// Frequency prewarping: dt_warp makes the bilinear map exact at the critical frequency w_c.
// @cite astrom1997 -- Astrom & Wittenmark, "Computer-Controlled Systems", 3rd ed., 1997, Sec. 3.5
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY> discretize(tustin_prewarp<Scalar> warp, const continuous_state_space<Scalar, NX, NU, NY>& sys, Scalar dt)
{
    const Scalar dt_warp = (Scalar{2} / warp.w_c) * std::tan(warp.w_c * dt / Scalar{2});
    return detail::tustin_bilinear(sys, dt_warp);
}

// @cite astrom1997 -- Astrom & Wittenmark, "Computer-Controlled Systems", 3rd ed., 1997, Sec. 3.5
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY> discretize(forward_euler, const continuous_state_space<Scalar, NX, NU, NY>& sys, Scalar dt)
{
    const Matrix<Scalar, NX, NX> identity = Matrix<Scalar, NX, NX>::Identity();

    Matrix<Scalar, NX, NX> Ad = identity + sys.A * dt;
    Matrix<Scalar, NX, NU> Bd = sys.B * dt;

    return {Ad, Bd, sys.C, sys.D};
}

// @cite astrom1997 -- Astrom & Wittenmark, "Computer-Controlled Systems", 3rd ed., 1997, Sec. 3.5
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY> discretize(backward_euler, const continuous_state_space<Scalar, NX, NU, NY>& sys, Scalar dt)
{
    const Matrix<Scalar, NX, NX> identity = Matrix<Scalar, NX, NX>::Identity();
    const Matrix<Scalar, NX, NX> backward_resolvent_inverse = (identity - sys.A * dt).inverse();

    Matrix<Scalar, NX, NX> Ad = backward_resolvent_inverse;
    Matrix<Scalar, NX, NU> Bd = backward_resolvent_inverse * sys.B * dt;
    Matrix<Scalar, NY, NX> Cd = sys.C * backward_resolvent_inverse;
    Matrix<Scalar, NY, NU> Dd = sys.D + Cd * sys.B * dt;

    return {Ad, Bd, Cd, Dd};
}

}

#endif
