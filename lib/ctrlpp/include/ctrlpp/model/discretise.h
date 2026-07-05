#ifndef HPP_GUARD_CTRLPP_MODEL_DISCRETISE_H
#define HPP_GUARD_CTRLPP_MODEL_DISCRETISE_H

/// @brief Continuous-to-discrete state-space conversion (ZOH, Tustin, Euler).
///
/// @cite franklin2015 -- Franklin et al., "Feedback Control of Dynamic Systems", 2015
/// @cite astrom1997 -- Astrom & Wittenmark, "Computer-Controlled Systems: Theory and Design", 3rd ed., 1997, Ch. 3

#include "ctrlpp/types.h"

#include "ctrlpp/model/state_space.h"

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

/// @brief Tustin (bilinear) discretisation with frequency prewarping.
///
/// Rescales the sample period so the bilinear map is exact at the chosen critical
/// frequency `w_c` (rad/s), trading pole-mapping accuracy elsewhere for exactness there.
template <typename Scalar>
struct tustin_prewarp
{
    Scalar w_c;
};

struct forward_euler
{
};

struct backward_euler
{
};

/// @brief Zero-order-hold discretisation via Van Loan's augmented matrix exponential method.
///
/// Forms the augmented matrix M = [[A*dt, B*dt], [0, 0]], computes exp(M), and extracts
/// Ad and Bd from the upper blocks. Cd = C, Dd = D.
///
/// @cite vanloan1978 -- Van Loan, "Computing Integrals Involving the Matrix Exponential",
///   IEEE Trans. Autom. Control 23(3):395-404, 1978
/// @cite astrom1997 -- Astrom & Wittenmark, "Computer-Controlled Systems", 3rd ed., 1997, Sec. 3.2
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY> discretise(zoh, const continuous_state_space<Scalar, NX, NU, NY>& sys, Scalar dt)
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);
    constexpr int aug = nx + nu;

    // Form augmented matrix
    Eigen::Matrix<Scalar, aug, aug> M;
    M.setZero();
    M.template block<nx, nx>(0, 0) = sys.A * dt;
    M.template block<nx, nu>(0, nx) = sys.B * dt;

    // Compute matrix exponential
    Eigen::Matrix<Scalar, aug, aug> expM = M.exp();

    // Extract discretised matrices
    Matrix<Scalar, NX, NX> Ad = expM.template block<nx, nx>(0, 0);
    Matrix<Scalar, NX, NU> Bd = expM.template block<nx, nu>(0, nx);

    return {Ad, Bd, sys.C, sys.D};
}

// Convenience wrapper matching the generic discretise<Method>(...) signature.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY> discretise(const continuous_state_space<Scalar, NX, NU, NY>& sys, Scalar dt, zoh = {})
{
    return discretise(zoh{}, sys, dt);
}

namespace detail
{

/// @brief Shared bilinear (Tustin) transform, parameterised on the (possibly prewarped) sample period.
///
/// Ad = (I - A*dt/2)^-1 (I + A*dt/2), Bd = (I - A*dt/2)^-1 B*dt, Cd = C (I - A*dt/2)^-1,
/// and the biproper feed-through correction Dd = D + C*Bd/2, which accounts for the
/// output coupling a bilinear map introduces even when the continuous system is strictly proper.
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

/// @brief Tustin (bilinear transform) discretisation.
///
/// @cite astrom1997 -- Astrom & Wittenmark, "Computer-Controlled Systems", 3rd ed., 1997, Sec. 3.5
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY> discretise(tustin, const continuous_state_space<Scalar, NX, NU, NY>& sys, Scalar dt)
{
    return detail::tustin_bilinear(sys, dt);
}

/// @brief Frequency-prewarped Tustin discretisation.
///
/// Replaces dt by dt_warp = (2/w_c) * tan(w_c*dt/2) before applying the bilinear map, so the
/// discrete and continuous frequency responses agree exactly at the critical frequency w_c.
///
/// @cite astrom1997 -- Astrom & Wittenmark, "Computer-Controlled Systems", 3rd ed., 1997, Sec. 3.5
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY> discretise(tustin_prewarp<Scalar> warp, const continuous_state_space<Scalar, NX, NU, NY>& sys, Scalar dt)
{
    const Scalar dt_warp = (Scalar{2} / warp.w_c) * std::tan(warp.w_c * dt / Scalar{2});
    return detail::tustin_bilinear(sys, dt_warp);
}

/// @brief Forward Euler discretisation.
///
/// Ad = I + A*dt, Bd = B*dt, Cd = C, Dd = D. First-order accurate; the discrete pole
/// z = 1 + s*dt is the first-order Taylor expansion of z = e^(s*dt) about s*dt = 0.
///
/// @cite astrom1997 -- Astrom & Wittenmark, "Computer-Controlled Systems", 3rd ed., 1997, Sec. 3.5
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY> discretise(forward_euler, const continuous_state_space<Scalar, NX, NU, NY>& sys, Scalar dt)
{
    const Matrix<Scalar, NX, NX> identity = Matrix<Scalar, NX, NX>::Identity();

    Matrix<Scalar, NX, NX> Ad = identity + sys.A * dt;
    Matrix<Scalar, NX, NU> Bd = sys.B * dt;

    return {Ad, Bd, sys.C, sys.D};
}

/// @brief Backward Euler discretisation.
///
/// Ad = (I - A*dt)^-1, Bd = (I - A*dt)^-1 B*dt, Cd = C (I - A*dt)^-1,
/// Dd = D + C (I - A*dt)^-1 B*dt. First-order accurate and unconditionally stable
/// for a stable continuous system, at the cost of a fixed-size matrix inversion.
///
/// @cite astrom1997 -- Astrom & Wittenmark, "Computer-Controlled Systems", 3rd ed., 1997, Sec. 3.5
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY> discretise(backward_euler, const continuous_state_space<Scalar, NX, NU, NY>& sys, Scalar dt)
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
