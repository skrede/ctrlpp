// These anchors run seeded Monte-Carlo ensembles of KF/EKF/UKF/MEKF against a
// ground-truth trajectory whose process and measurement noise are sampled
// from exactly the covariances (Q, R) each filter is configured with, then
// check that the filter's own reported covariance is a statistically honest
// description of its actual error. The standard diagnostics for this are the
// Normalized Estimation Error Squared (NEES, using ground truth) and the
// Normalized Innovation Squared (NIS, using only the filter's own residual
// and innovation covariance); both are chi-square distributed under a
// correctly tuned, correctly implemented filter (Bar-Shalom, Li & Kirubarajan,
// "Estimation with Applications to Tracking and Navigation", 2001, Sec. 5.4).
//
// The acceptance band for the average of M independent chi-square(k) samples
// is derived analytically: the sum of M iid chi-square(k) variables is itself
// chi-square(M*k) distributed, so a two-sided (1-alpha) confidence interval
// for the average is [chi2_quantile(alpha/2, M*k), chi2_quantile(1-alpha/2,
// M*k)] / M. The chi-square quantile itself is evaluated via the Wilson-Hilferty
// cube-root normal approximation (Wilson & Hilferty, "The distribution of
// chi-square", PNAS 17, 1931), which is highly accurate for the large degrees
// of freedom (M*k, with M in the thousands) used here; the two-sided 95%
// level uses the standard normal 97.5th percentile z=1.959964 (Abramowitz &
// Stegun, "Handbook of Mathematical Functions", 1964, Table 26.1).
//
// KF, EKF, UKF, and MEKF are all correct now: every section is untagged and
// must pass. The UKF's posterior covariance update is the algebraically
// complete P - K*S*K^T reduction with no spurious K*R*K^T inflation, and its
// sigma points are built from a square root that reconstructs the covariance
// exactly. The MEKF's error-state transition propagates the attitude sub-block
// with the transpose of the incremental rotation, matching its right-error
// multiplicative-correction convention (invisible on isotropic covariance,
// exposed here once the attitude covariance is anisotropic). The MEKF section
// samples its initial true attitude error from P0 (Sec. 5.4) so the single
// unobservable rotational DOF of the vector measurement does not bias the NEES
// low; with the corrected transition its NEES lies inside the band.

#include "ctrlpp/lie/so3.h"
#include "ctrlpp/estimation/ekf.h"
#include "ctrlpp/estimation/ukf.h"
#include "ctrlpp/estimation/mekf.h"
#include "ctrlpp/estimation/kalman.h"
#include "ctrlpp/model/state_space.h"
#include "ctrlpp/estimation/sigma_points/merwe_sigma_points.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <random>
#include <cstddef>

using namespace ctrlpp;

namespace
{

/// Wilson-Hilferty cube-root normal approximation to the chi-square quantile
/// function, accurate to a small fraction of a percent for the degrees of
/// freedom used in this file (all in the thousands).
/// @cite wilson1931 -- Wilson & Hilferty, "The distribution of chi-square", PNAS 17, 684-688, 1931
double chi_square_quantile(double dof, double z)
{
    double term = 1.0 - 2.0 / (9.0 * dof) + z * std::sqrt(2.0 / (9.0 * dof));
    return dof * term * term * term;
}

/// Standard normal 97.5th percentile: Phi^-1(0.975).
/// @cite abramowitz1964 -- Abramowitz & Stegun, "Handbook of Mathematical Functions", 1964, Table 26.1
constexpr double z_975 = 1.9599639845400545;

struct chi_square_band
{
    double lower;
    double upper;
};

/// Two-sided 95% confidence band for the average of `num_samples` iid
/// chi-square(dim) variables, derived from their sum being chi-square(dim *
/// num_samples) distributed.
chi_square_band average_chi_square_band(std::size_t dim, std::size_t num_samples)
{
    double dof = static_cast<double>(dim) * static_cast<double>(num_samples);
    double n = static_cast<double>(num_samples);
    return {chi_square_quantile(dof, -z_975) / n, chi_square_quantile(dof, z_975) / n};
}

/// Draw a zero-mean Gaussian vector with covariance `cov` via its Cholesky factor.
template <std::size_t N>
Vector<double, N> sample_gaussian(const Matrix<double, N, N>& cov, std::mt19937& gen)
{
    std::normal_distribution<double> nd(0.0, 1.0);
    Vector<double, N> z;
    for(std::size_t i = 0; i < N; ++i)
        z(static_cast<Eigen::Index>(i)) = nd(gen);
    Eigen::LLT<Matrix<double, N, N>> llt(cov);
    return (llt.matrixL() * z).eval();
}

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr std::size_t NY = 1;
constexpr double dt = 0.1;
constexpr std::size_t T = 20;
constexpr std::size_t M = 2000;

struct linear_dynamics
{
    auto operator()(const Vector<double, NX>& x, const Vector<double, NU>& u) const -> Vector<double, NX>
    {
        Vector<double, NX> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    }
};

struct position_measurement
{
    auto operator()(const Vector<double, NX>& x) const -> Vector<double, NY>
    {
        return (Vector<double, NY>() << x(0)).finished();
    }
};

/// Mildly nonlinear dynamics: a bounded pendulum-like restoring term added to
/// the linear kinematics above. Its curvature is bounded by k_nl regardless of
/// state magnitude, so its linearization error stays a small fraction of the
/// per-step process noise over the horizon used here, keeping the EKF close
/// to its own linear-Gaussian ideal.
struct mild_nonlinear_dynamics
{
    static constexpr double k_nl = 0.2;

    auto operator()(const Vector<double, NX>& x, const Vector<double, NU>& u) const -> Vector<double, NX>
    {
        Vector<double, NX> x_next;
        x_next(0) = x(0) + dt * x(1);
        x_next(1) = x(1) + dt * u(0) - dt * k_nl * std::sin(x(0));
        return x_next;
    }

    auto jacobian_x(const Vector<double, NX>& x, const Vector<double, NU>& /*u*/) const -> Matrix<double, NX, NX>
    {
        Matrix<double, NX, NX> F;
        F << 1.0, dt, -dt * k_nl * std::cos(x(0)), 1.0;
        return F;
    }

    auto jacobian_u(const Vector<double, NX>& /*x*/, const Vector<double, NU>& /*u*/) const -> Matrix<double, NX, NU>
    {
        Matrix<double, NX, NU> G;
        G << 0.0, dt;
        return G;
    }
};

using UkfType = ukf<double, NX, NU, NY, linear_dynamics, position_measurement>;
using EkfType = ekf<double, NX, NU, NY, mild_nonlinear_dynamics, position_measurement>;

}

TEST_CASE("KF NEES/NIS Monte-Carlo average lies within the chi-square consistency band", "[kalman][anchor]")
{
    discrete_state_space<double, NX, NU, NY> sys;
    sys.A << 1.0, dt, 0.0, 1.0;
    sys.B << 0.5 * dt * dt, dt;
    sys.C << 1.0, 0.0;
    sys.D = Matrix<double, NY, NU>::Zero();

    kalman_config<double, NX, NU, NY> cfg;
    cfg.Q = Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = Matrix<double, NY, NY>::Identity() * 0.1;
    cfg.P0 = Matrix<double, NX, NX>::Identity() * 10.0;

    std::mt19937 gen(20260705);
    Vector<double, NU> u = Vector<double, NU>::Zero();

    double nees_sum = 0.0;
    double nis_sum = 0.0;
    for(std::size_t m = 0; m < M; ++m)
    {
        kalman_filter<double, NX, NU, NY> filt(sys, cfg);
        Vector<double, NX> x_true = Vector<double, NX>::Zero();

        for(std::size_t t = 0; t < T; ++t)
        {
            Vector<double, NX> w = sample_gaussian<NX>(cfg.Q, gen);
            x_true = (sys.A * x_true + sys.B * u + w).eval();

            filt.predict(u);

            Vector<double, NY> v = sample_gaussian<NY>(cfg.R, gen);
            Vector<double, NY> z = (sys.C * x_true + v).eval();
            filt.update(z);
        }

        Vector<double, NX> e = (x_true - filt.state()).eval();
        Vector<double, NX> Pinv_e = filt.covariance().ldlt().solve(e);
        nees_sum += (e.transpose() * Pinv_e)(0, 0);
        nis_sum += filt.nis(); // filter's Normalized Innovation Squared: innovation' * S^-1 * innovation
    }

    double nees_avg = nees_sum / static_cast<double>(M);
    double nis_avg = nis_sum / static_cast<double>(M);
    auto nees_band = average_chi_square_band(NX, M);
    auto nis_band = average_chi_square_band(NY, M);

    CAPTURE(nees_avg, nees_band.lower, nees_band.upper);
    REQUIRE(nees_avg >= nees_band.lower);
    REQUIRE(nees_avg <= nees_band.upper);

    CAPTURE(nis_avg, nis_band.lower, nis_band.upper);
    REQUIRE(nis_avg >= nis_band.lower);
    REQUIRE(nis_avg <= nis_band.upper);
}

TEST_CASE("EKF NEES/NIS Monte-Carlo average lies within the chi-square consistency band on a mild nonlinearity", "[ekf][anchor]")
{
    ekf_config<double, NX, NU, NY> cfg;
    cfg.Q = Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = Matrix<double, NY, NY>::Identity() * 0.1;
    cfg.P0 = Matrix<double, NX, NX>::Identity() * 10.0;

    std::mt19937 gen(20260705);
    Vector<double, NU> u = Vector<double, NU>::Zero();
    mild_nonlinear_dynamics dyn;
    position_measurement meas;

    double nees_sum = 0.0;
    double nis_sum = 0.0;
    for(std::size_t m = 0; m < M; ++m)
    {
        EkfType filt(dyn, meas, cfg);
        Vector<double, NX> x_true = Vector<double, NX>::Zero();

        for(std::size_t t = 0; t < T; ++t)
        {
            Vector<double, NX> w = sample_gaussian<NX>(cfg.Q, gen);
            x_true = (dyn(x_true, u) + w).eval();

            filt.predict(u);

            Vector<double, NY> v = sample_gaussian<NY>(cfg.R, gen);
            Vector<double, NY> z = (meas(x_true) + v).eval();
            filt.update(z);
        }

        Vector<double, NX> e = (x_true - filt.state()).eval();
        Vector<double, NX> Pinv_e = filt.covariance().ldlt().solve(e);
        nees_sum += (e.transpose() * Pinv_e)(0, 0);
        nis_sum += filt.nis();
    }

    double nees_avg = nees_sum / static_cast<double>(M);
    double nis_avg = nis_sum / static_cast<double>(M);
    auto nees_band = average_chi_square_band(NX, M);
    auto nis_band = average_chi_square_band(NY, M);

    CAPTURE(nees_avg, nees_band.lower, nees_band.upper);
    REQUIRE(nees_avg >= nees_band.lower);
    REQUIRE(nees_avg <= nees_band.upper);

    CAPTURE(nis_avg, nis_band.lower, nis_band.upper);
    REQUIRE(nis_avg >= nis_band.lower);
    REQUIRE(nis_avg <= nis_band.upper);
}

TEST_CASE("UKF NEES/NIS Monte-Carlo average lies within the chi-square consistency band", "[ukf][anchor]")
{
    // With the posterior covariance update reduced to P - K*S*K^T and the
    // sigma-point square root reconstructing the covariance exactly, the
    // reported covariance tracks the true error statistics, so both the NEES
    // and NIS averages fall inside the chi-square consistency band.
    merwe_options<double> strategy_opts;
    strategy_opts.alpha = 1.0;
    strategy_opts.beta = 0.0;
    strategy_opts.kappa = 3.0 - static_cast<double>(NX);

    ukf_config<double, NX, NU, NY> cfg;
    cfg.Q = Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = Matrix<double, NY, NY>::Identity() * 0.1;
    cfg.P0 = Matrix<double, NX, NX>::Identity() * 10.0;

    Matrix<double, NX, NX> A;
    A << 1.0, dt, 0.0, 1.0;
    Matrix<double, NX, NU> B;
    B << 0.5 * dt * dt, dt;
    Matrix<double, NY, NX> C;
    C << 1.0, 0.0;

    std::mt19937 gen(20260705);
    Vector<double, NU> u = Vector<double, NU>::Zero();

    double nees_sum = 0.0;
    double nis_sum = 0.0;
    for(std::size_t m = 0; m < M; ++m)
    {
        UkfType filt(linear_dynamics{}, position_measurement{}, cfg, strategy_opts);
        Vector<double, NX> x_true = Vector<double, NX>::Zero();
        Matrix<double, NX, NX> P_prior = cfg.P0;

        for(std::size_t t = 0; t < T; ++t)
        {
            Vector<double, NX> w = sample_gaussian<NX>(cfg.Q, gen);
            x_true = (A * x_true + B * u + w).eval();

            filt.predict(u);
            P_prior = filt.covariance();

            Vector<double, NY> v = sample_gaussian<NY>(cfg.R, gen);
            Vector<double, NY> z = (C * x_true + v).eval();
            filt.update(z);
        }

        Vector<double, NX> e = (x_true - filt.state()).eval();
        Vector<double, NX> Pinv_e = filt.covariance().ldlt().solve(e);
        nees_sum += (e.transpose() * Pinv_e)(0, 0);

        // Reconstruct the innovation covariance S = C*P_prior*C^T + R ourselves
        // from the known linear measurement map and the filter's own reported
        // prior covariance, matching how the KF/EKF's own internal S is defined.
        // This independent reconstruction cross-checks the sigma-point path
        // rather than reusing the filter's nis() accessor.
        Matrix<double, NY, NY> S = (C * P_prior * C.transpose() + cfg.R).eval();
        Vector<double, NY> innovation = filt.innovation();
        Vector<double, NY> Sinv_innovation = S.ldlt().solve(innovation);
        nis_sum += (innovation.transpose() * Sinv_innovation)(0, 0);
    }

    double nees_avg = nees_sum / static_cast<double>(M);
    double nis_avg = nis_sum / static_cast<double>(M);
    auto nees_band = average_chi_square_band(NX, M);
    auto nis_band = average_chi_square_band(NY, M);

    CAPTURE(nees_avg, nees_band.lower, nees_band.upper);
    REQUIRE(nees_avg >= nees_band.lower);
    REQUIRE(nees_avg <= nees_band.upper);

    CAPTURE(nis_avg, nis_band.lower, nis_band.upper);
    REQUIRE(nis_avg >= nis_band.lower);
    REQUIRE(nis_avg <= nis_band.upper);
}

TEST_CASE("MEKF attitude NEES Monte-Carlo average lies within the chi-square consistency band on an anisotropic P", "[estimation][anchor]")
{
    // The error-state transition propagates the attitude sub-block with the
    // transpose of the incremental rotation, matching the filter's right-error
    // multiplicative-correction convention (apply_multiplicative_correction:
    // q_ = q_ * exp(delta_att)). The effect is invisible on isotropic covariance
    // (a rotation commutes with a scalar multiple of identity), so an anisotropic
    // initial attitude covariance is used here to expose it; with the corrected
    // transition the attitude NEES lands inside the chi-square consistency band.
    constexpr std::size_t NB = 3; // mekf's predict_impl requires a >=3-dim bias (b_.head<3>())
    constexpr std::size_t NY_MEKF = 3;

    struct vector_observation_measurement
    {
        Vector<double, 3> r_w;
        auto operator()(const Eigen::Quaternion<double>& q, const Vector<double, NB>& /*b*/) const -> Vector<double, NY_MEKF>
        {
            return q.conjugate() * r_w;
        }
    };

    Vector<double, 3> r_w;
    r_w << 0.0, 0.0, 1.0;

    mekf_config<double, NB, NY_MEKF> cfg;
    cfg.P0.setZero();
    // Anisotropic (distinct per-axis variances) attitude block: the isotropic
    // blind spot documented in the anisotropic-P covariance-axis anchor.
    cfg.P0.template block<3, 3>(0, 0) = Vector<double, 3>{0.09, 0.01, 0.0025}.asDiagonal();
    cfg.P0.template block<NB, NB>(3, 3) = Matrix<double, NB, NB>::Identity() * 1e-6;
    cfg.Q.setZero();
    cfg.Q.template block<3, 3>(0, 0) = Matrix<double, 3, 3>::Identity() * 1e-4;
    cfg.Q.template block<NB, NB>(3, 3) = Matrix<double, NB, NB>::Identity() * 1e-10;
    cfg.R = Matrix<double, NY_MEKF, NY_MEKF>::Identity() * 0.01;
    cfg.dt = dt;

    Vector<double, 3> omega_true;
    omega_true << 0.2, -0.1, 0.15;

    std::mt19937 gen(20260705);
    vector_observation_measurement meas{r_w};

    double nees_sum = 0.0;
    for(std::size_t m = 0; m < M; ++m)
    {
        mekf<double, NB, NY_MEKF, vector_observation_measurement> filt(meas, cfg);
        // Sample the initial true attitude error from P0 so the NEES is
        // consistent from t=0 (Bar-Shalom, Li & Kirubarajan 2001, Sec. 5.4).
        // The single-vector measurement leaves one rotational DOF unobservable,
        // so a zero initial error against a nonzero P0 would never wash out and
        // would bias the NEES low. The filter estimate starts at q0 = Identity,
        // so under the right-error convention q_true = q_est * exp(dtheta0).
        Vector<double, 3> dtheta0 = sample_gaussian<3>(cfg.P0.template block<3, 3>(0, 0), gen);
        Eigen::Quaternion<double> q_true = so3::exp(dtheta0);

        for(std::size_t t = 0; t < T; ++t)
        {
            Vector<double, 3> eta = sample_gaussian<3>(cfg.Q.template block<3, 3>(0, 0), gen);
            Vector<double, 3> omega_dt = (omega_true * dt + eta).eval();
            q_true = (q_true * so3::exp(omega_dt)).normalized();

            filt.predict(omega_true, dt);

            Vector<double, NY_MEKF> v = sample_gaussian<NY_MEKF>(cfg.R, gen);
            Vector<double, NY_MEKF> z = (q_true.conjugate() * r_w + v).eval();
            filt.update(z);
        }

        Eigen::Quaternion<double> q_hat = filt.attitude();
        Vector<double, 3> e_att = so3::log((q_hat.conjugate() * q_true).normalized());
        Matrix<double, 3, 3> P_att = filt.covariance().template block<3, 3>(0, 0);
        Vector<double, 3> Pinv_e = P_att.ldlt().solve(e_att);
        nees_sum += (e_att.transpose() * Pinv_e)(0, 0);
    }

    double nees_avg = nees_sum / static_cast<double>(M);
    auto nees_band = average_chi_square_band(3, M);

    CAPTURE(nees_avg, nees_band.lower, nees_band.upper);
    REQUIRE(nees_avg >= nees_band.lower);
    REQUIRE(nees_avg <= nees_band.upper);
}
