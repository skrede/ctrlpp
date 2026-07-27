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
// guarantees by definition. The generator's matrix square root is now the
// unpivoted Cholesky factor, which satisfies S*S^T = P exactly with no
// permutation to track, so the reconstruction squares back to the original P
// even for anisotropic matrices that would force a pivoted factorization to
// permute. A companion case checks that the SO(3) manifold sigma points,
// which delegate to the same square root, inherit the fix on the same
// anisotropic tangent-space covariance.
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
// The sigma-point cases now pass as active tests. The MEKF case's transpose
// defect is not yet corrected, so it stays tagged with [!shouldfail] until
// that transition matrix is fixed, at which point its tag is removed too.

#include "ctrlpp/lie/so3.h"
#include "ctrlpp/estimation/mekf.h"
#include "ctrlpp/estimation/sigma_points/merwe_sigma_points.h"
#include "ctrlpp/estimation/sigma_points/so3_sigma_points.h"

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

TEST_CASE("unscented sigma points reconstruct a non-identity-pivot covariance", "[estimation][anchor]")
{
    constexpr std::size_t SP_NX = 3;

    merwe_options<double> opts;
    opts.alpha = 1.0;
    opts.beta = 0.0;
    opts.kappa = 3.0 - static_cast<double>(SP_NX);

    auto strategy_result = merwe_sigma_points<double, SP_NX>::try_create(opts);
    REQUIRE(strategy_result.has_value());
    auto const& strategy = *strategy_result;

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

TEST_CASE("SO(3) manifold sigma points reconstruct a non-identity-pivot tangent covariance", "[estimation][anchor]")
{
    merwe_options<double> opts;
    opts.alpha = 1.0;
    opts.beta = 0.0;
    opts.kappa = 3.0 - 3.0;

    auto strategy_result = so3_merwe_sigma_points<double>::try_create(opts);
    REQUIRE(strategy_result.has_value());
    auto const& strategy = *strategy_result;

    // The same anisotropic structure, scaled so that every sigma-point offset
    // stays well inside the SO(3) exponential's injectivity radius (norm below
    // pi); at that scale the logarithm of each composed point recovers its
    // tangent offset exactly, isolating the square-root reconstruction that the
    // manifold generator inherits from the flat one.
    const Matrix<double, 3, 3> P = (anisotropic_covariance_3x3() * 0.01).eval();
    const Eigen::Quaternion<double> q_mean = Eigen::Quaternion<double>::Identity();

    auto sigma = strategy.generate(q_mean, P);

    // The manifold sigma points sit at q_mean composed with the exponential of
    // the tangent-space offsets; with q_mean the identity the logarithm of each
    // point recovers its tangent offset, so the weighted covariance of those
    // logarithms reconstructs the tangent covariance P.
    Matrix<double, 3, 3> reconstructed = Matrix<double, 3, 3>::Zero();
    for(std::size_t i = 0; i < sigma.points.size(); ++i)
    {
        Vector<double, 3> tangent = so3::log(sigma.points[i]);
        reconstructed += sigma.Wc[i] * tangent * tangent.transpose();
    }

    // The reconstruction sums a transcendental exp/log round trip over all
    // 2*NX+1 sigma points, so the accumulated rounding scales with the number
    // of points and the covariance magnitude.
    const double eps = std::numeric_limits<double>::epsilon();
    const double tol = static_cast<double>(so3_merwe_sigma_points<double>::num_points) * eps * (1.0 + P.norm());

    for(std::size_t i = 0; i < 3; ++i)
    {
        for(std::size_t j = 0; j < 3; ++j)
        {
            CAPTURE(i, j);
            REQUIRE_THAT(reconstructed(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)),
                WithinAbs(P(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)), tol));
        }
    }
}

TEST_CASE("MEKF error-state transition matches the analytic covariance transform on an anisotropic P", "[estimation][anchor]")
{
    mekf_config<double, NB, NY> cfg;
    cfg.P0.setZero();
    cfg.P0.template block<3, 3>(0, 0) = anisotropic_covariance_3x3();
    cfg.P0.template block<nb, nb>(3, 3) = Matrix<double, NB, NB>::Identity() * 1e-4;
    cfg.Q.setZero();
    cfg.dt = 0.1;

    auto filter_result = mekf<double, NB, NY, trivial_mekf_measurement>::create(trivial_mekf_measurement{}, cfg);
    REQUIRE(filter_result.has_value());
    auto& filter = *filter_result;

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
