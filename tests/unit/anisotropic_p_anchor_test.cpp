// These anchors probe the one blind spot that isotropic (scalar-multiple-
// of-identity) covariance tests cannot see: any bug that only shows up when
// a covariance matrix has distinct eigenvalues and off-diagonal
// correlation. A rotation R and a coordinate-aligned identity or uniform
// covariance commute trivially (R * (c*I) * R^T = c*I for any orthogonal
// R), so a transpose or permutation error hidden inside a rotation-like
// operator produces zero error on isotropic P and only appears once P is
// genuinely anisotropic.
//
// The first case constructs a strongly anisotropic, non-identity-pivot
// covariance and checks that the unscented sigma point set's weighted
// covariance reconstructs it exactly, as the sigma point construction
// guarantees by definition. It currently fails because the sigma point
// generator's matrix square root is read off an LDLT factorization while
// discarding that factorization's row/column pivot permutation: the
// factorization satisfies P(permuted) = L*D*L^T, not P = L*D*L^T, so the
// naive L*sqrt(D) read-off only squares back to the original P when no
// pivoting occurred, which anisotropic matrices routinely require.
//
// The second case checks the MEKF's one-step error-state covariance
// transition against an independently computed analytic transform on the
// same anisotropic attitude block. It currently fails because the
// transition uses the incremental rotation where its own multiplicative
// correction convention requires the incremental rotation's transpose;
// since a rotation matrix is orthogonal, this transpose error vanishes
// exactly on isotropic P and only surfaces once the attitude covariance is
// anisotropic.
//
// Both REQUIRE blocks below currently fail against pre-fix code;
// [!shouldfail] reports each case as passing until its defect is
// corrected, at which point the tag on that case must be removed.

#include "ctrlpp/lie/so3.h"
#include "ctrlpp/estimation/mekf.h"
#include "ctrlpp/estimation/sigma_points/merwe_sigma_points.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <limits>
#include <cstddef>

using namespace ctrlpp;
using Catch::Matchers::WithinAbs;

namespace
{

/// @brief A strongly anisotropic, non-identity-pivot 3x3 covariance
/// (eigenvalues 100, 4, 0.01 under a generic rotation), constructed so
/// that its LDLT factorization requires a genuine pivot permutation.
auto anisotropic_covariance_3x3() -> Matrix<double, 3, 3>
{
    Matrix<double, 3, 3> P;
    P << 4.6030568403028411, 8.2110257815469581, -13.971872104544538,
        8.2110257815469581, 37.551437765319832, -46.007225698288167,
        -13.971872104544538, -46.007225698288167, 61.855505394377325;
    return P;
}

constexpr std::size_t NB = 3;
constexpr std::size_t NY = 1;
constexpr std::size_t NE = 3 + NB;
constexpr int nb = static_cast<int>(NB);

struct trivial_mekf_measurement
{
    auto operator()(const Eigen::Quaternion<double>&, const Vector<double, NB>&) const -> Vector<double, NY>
    {
        return Vector<double, NY>::Zero();
    }
};

} // namespace

TEST_CASE("unscented sigma points reconstruct a non-identity-pivot covariance", "[estimation][anchor][!shouldfail]")
{
    constexpr std::size_t SP_NX = 3;

    merwe_options<double> opts;
    opts.alpha = 1.0;
    opts.beta = 0.0;
    opts.kappa = 3.0 - static_cast<double>(SP_NX);

    merwe_sigma_points<double, SP_NX> strategy(opts);

    const Matrix<double, SP_NX, SP_NX> P = anisotropic_covariance_3x3();
    const Vector<double, SP_NX> x0 = Vector<double, SP_NX>::Zero();

    auto sigma = strategy.generate(x0, P);

    Matrix<double, SP_NX, SP_NX> reconstructed = Matrix<double, SP_NX, SP_NX>::Zero();
    for(std::size_t i = 0; i < sigma.points.size(); ++i)
    {
        auto diff = (sigma.points[i] - x0).eval();
        reconstructed += sigma.Wc[i] * diff * diff.transpose();
    }

    // The weighted sigma point covariance reproduces P exactly by
    // construction (Wan & van der Merwe 2001, Eq. 15/18); the only
    // admissible error is floating-point rounding, which scales with
    // problem size and the covariance's own magnitude.
    const double eps = std::numeric_limits<double>::epsilon();
    const double tol = static_cast<double>(SP_NX) * eps * (1.0 + P.norm());

    for(std::size_t i = 0; i < SP_NX; ++i)
    {
        for(std::size_t j = 0; j < SP_NX; ++j)
        {
            CAPTURE(i, j);
            REQUIRE_THAT(reconstructed(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)),
                WithinAbs(P(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)), tol));
        }
    }
}

TEST_CASE("MEKF error-state transition matches the analytic covariance transform on an anisotropic P", "[estimation][anchor][!shouldfail]")
{
    mekf_config<double, NB, NY> cfg;
    cfg.P0.setZero();
    cfg.P0.template block<3, 3>(0, 0) = anisotropic_covariance_3x3();
    cfg.P0.template block<nb, nb>(3, 3) = Matrix<double, NB, NB>::Identity() * 1e-4;
    cfg.Q.setZero();
    cfg.dt = 0.1;

    mekf<double, NB, NY, trivial_mekf_measurement> filter(trivial_mekf_measurement{}, cfg);

    Vector<double, 3> omega;
    omega << 0.3, -0.2, 0.5;

    filter.predict(omega, cfg.dt);

    // Independent analytic oracle: for a right-error (body-frame)
    // multiplicative attitude error, the correct error-state transition
    // propagates the attitude sub-block with the transpose of the
    // incremental rotation (Markley, "Attitude Error Representations for
    // Kalman Filtering", 2003; Sola et al., "A micro Lie theory for state
    // estimation in robotics", 2018). This is identical to the filter's
    // own transition matrix in every respect except that one transpose.
    Vector<double, 3> omega_dt = (omega * cfg.dt).eval();
    Matrix<double, 3, 3> C = so3::exp(omega_dt).toRotationMatrix();

    Matrix<double, NE, NE> F_correct = Matrix<double, NE, NE>::Identity();
    F_correct.template block<3, 3>(0, 0) = C.transpose();
    F_correct.template block<3, nb>(0, 3) = -Matrix<double, 3, NB>::Identity() * cfg.dt;

    Matrix<double, NE, NE> P_correct = detail::symmetrize((F_correct * cfg.P0 * F_correct.transpose()).eval());

    const double eps = std::numeric_limits<double>::epsilon();
    const double tol = static_cast<double>(NE) * eps * (1.0 + P_correct.norm());

    for(std::size_t i = 0; i < NE; ++i)
    {
        for(std::size_t j = 0; j < NE; ++j)
        {
            CAPTURE(i, j);
            REQUIRE_THAT(filter.covariance()(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)),
                WithinAbs(P_correct(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)), tol));
        }
    }
}
