#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_SMOOTHING_SPLINE_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_SMOOTHING_SPLINE_H

/// @brief Smoothing spline with configurable mu tradeoff parameter.
///
/// Constructs a C2-continuous smoothing spline that approximates noisy waypoint
/// data. The mu parameter controls the tradeoff between data fidelity and
/// smoothness: mu=1 yields exact interpolation (passes through all waypoints),
/// while mu near 0 maximizes smoothness at the cost of data fit.
///
/// The regularization weight lambda = 2*(1-mu)/(3*mu) maps mu in (0,1] to
/// lambda in [0, inf). The system (R + lambda * Q^T * Q) * d = Q^T * q is
/// solved for second derivatives at interior knots, then smoothed positions
/// and cubic polynomial coefficients are computed per span.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
/// Automatic Machines and Robots", 2009, Sec. 4.4.5

#include "ctrlpp/trajectory/trajectory_segment.h"
#include "ctrlpp/trajectory/trajectory_types.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <vector>

namespace ctrlpp
{

/// @brief Smoothing spline approximation with mu tradeoff parameter.
///
/// mu=1.0 degenerates to interpolation (passes through all waypoints).
/// mu near 0 biases toward smoothness (deviates from data).
/// Natural-like endpoint conditions: d_0 = d_n = 0.
///
/// @cite biagiotti2009 -- Sec. 4.4.5
template <typename Scalar>
class smoothing_spline
{
  public:
    struct config
    {
        std::vector<Scalar> times;       ///< Knot times t_0 ... t_n (n+1 entries)
        std::vector<Scalar> positions;   ///< Waypoint positions q_0 ... q_n (n+1 entries)
        Scalar mu{Scalar{0.5}};          ///< Tradeoff: 1.0 = interpolation, near 0 = max smoothness
    };

    /// @brief Construct smoothing spline from waypoints and mu parameter.
    ///
    /// Solves the regularized system for second derivatives, computes smoothed
    /// positions, then derives cubic polynomial coefficients per span.
    ///
    /// @cite biagiotti2009 -- Sec. 4.4.5
    explicit smoothing_spline(config const& cfg)
        : times_{cfg.times}
    {
        auto const n_pts = times_.size();
        assert(n_pts >= 2);
        assert(cfg.positions.size() == n_pts);

        auto const n = n_pts - 1; // number of spans

        // Clamp mu to [epsilon, 1.0] to prevent degenerate lambda (Pitfall 2)
        auto const eps = Scalar{1e-10};
        auto const mu = std::clamp(cfg.mu, eps, Scalar{1});

        // Regularization weight: lambda = 2*(1-mu)/(3*mu)
        // @cite biagiotti2009 -- Sec. 4.4.5
        auto const lambda = Scalar{2} * (Scalar{1} - mu) / (Scalar{3} * mu);

        // Compute span durations h_i
        std::vector<Scalar> h(n);
        for (std::size_t i = 0; i < n; ++i) {
            h[i] = times_[i + 1] - times_[i];
            assert(h[i] > Scalar{0});
        }

        // For n_pts == 2, no interior knots -- degenerate to linear segment
        if (n_pts == 2) {
            build_linear(cfg.positions, h);
            return;
        }

        // Build and solve the system for interior second derivatives d_1..d_{n-1}
        // Natural-like endpoints: d_0 = d_n = 0
        // System: (R + lambda * Q^T * Q) * d = Q^T * q
        // where R is (n-1)x(n-1) tridiagonal, Q is (n+1)x(n-1) second-difference
        // @cite biagiotti2009 -- Sec. 4.4.5
        auto const m = n - 1; // interior knot count

        solve_and_build_coeffs(cfg.positions, h, lambda, n_pts, n, m);
    }

    /// @brief Evaluate smoothing spline at time t, clamped to [t_0, t_n].
    ///
    /// Uses binary search to find the active span, then Horner evaluation
    /// for position, velocity, and acceleration.
    ///
    /// @cite biagiotti2009 -- Sec. 4.4.5
    auto evaluate(Scalar t) const -> trajectory_point<Scalar, 1>
    {
        auto const tc = std::clamp(t, times_.front(), times_.back());
        auto const span = find_span(tc);

        auto const s = tc - times_[span];
        auto const& c = coeffs_[span];

        // Horner evaluation: q = a + s*(b + s*(c + s*d))
        auto const q = c[0] + s * (c[1] + s * (c[2] + s * c[3]));
        auto const dq = c[1] + s * (Scalar{2} * c[2] + Scalar{3} * c[3] * s);
        auto const ddq = Scalar{2} * c[2] + Scalar{6} * c[3] * s;

        return {
            .position = Vector<Scalar, 1>{q},
            .velocity = Vector<Scalar, 1>{dq},
            .acceleration = Vector<Scalar, 1>{ddq},
        };
    }

    /// @brief Total spline duration: t_n - t_0.
    auto duration() const -> Scalar { return times_.back() - times_.front(); }

  private:
    /// @brief Build linear segment for 2-point case.
    void build_linear(std::vector<Scalar> const& pos,
                      std::vector<Scalar> const& h)
    {
        auto const slope = (pos[1] - pos[0]) / h[0];
        coeffs_.resize(1);
        coeffs_[0][0] = pos[0];
        coeffs_[0][1] = slope;
        coeffs_[0][2] = Scalar{0};
        coeffs_[0][3] = Scalar{0};
    }

    /// @brief Solve regularized system and compute cubic coefficients.
    ///
    /// Builds (R + lambda * Q^T * Q) for the interior second derivatives,
    /// solves via Eigen dense QR, then computes smoothed positions and
    /// cubic polynomial coefficients per span.
    ///
    /// @cite biagiotti2009 -- Sec. 4.4.5
    void solve_and_build_coeffs(std::vector<Scalar> const& pos,
                                std::vector<Scalar> const& h,
                                Scalar lambda,
                                std::size_t n_pts,
                                std::size_t n,
                                std::size_t m)
    {
        using MatX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        using VecX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

        auto const mi = static_cast<Eigen::Index>(m);
        auto const ni = static_cast<Eigen::Index>(n_pts);

        // Build R: (m x m) tridiagonal matrix
        // R[i][i] = (h_{i} + h_{i+1}) / 3  (using 1-based interior indexing -> h[i], h[i+1])
        // R[i][i+1] = R[i+1][i] = h_{i+1} / 6
        // @cite biagiotti2009 -- Sec. 4.4.5
        MatX R = MatX::Zero(mi, mi);
        for (Eigen::Index i = 0; i < mi; ++i) {
            auto const ii = static_cast<std::size_t>(i);
            R(i, i) = (h[ii] + h[ii + 1]) / Scalar{3};
            if (i + 1 < mi) {
                R(i, i + 1) = h[ii + 1] / Scalar{6};
                R(i + 1, i) = h[ii + 1] / Scalar{6};
            }
        }

        // Build Q: (n_pts x m) second-difference matrix
        // Q[i][j] relates position i to interior second derivative j+1
        // For interior knot j (0-based), affecting positions j, j+1, j+2:
        //   Q[j][j]   = 1/h[j]
        //   Q[j+1][j] = -(1/h[j] + 1/h[j+1])
        //   Q[j+2][j] = 1/h[j+1]
        // @cite biagiotti2009 -- Sec. 4.4.5
        MatX Q = MatX::Zero(ni, mi);
        for (Eigen::Index j = 0; j < mi; ++j) {
            auto const jj = static_cast<std::size_t>(j);
            Q(j, j) = Scalar{1} / h[jj];
            Q(j + 1, j) = -(Scalar{1} / h[jj] + Scalar{1} / h[jj + 1]);
            Q(j + 2, j) = Scalar{1} / h[jj + 1];
        }

        // Build system matrix A = R + lambda * Q^T * Q
        MatX const QtQ = Q.transpose() * Q;
        MatX const A = R + lambda * QtQ;

        // Build RHS: Q^T * q
        VecX q(ni);
        for (Eigen::Index i = 0; i < ni; ++i) {
            q(i) = pos[static_cast<std::size_t>(i)];
        }
        VecX const rhs = Q.transpose() * q;

        // Solve for interior second derivatives d_1..d_{m}
        VecX const d_int = A.colPivHouseholderQr().solve(rhs);

        // Full second derivative vector (with d_0 = d_n = 0)
        std::vector<Scalar> d(n_pts, Scalar{0});
        for (Eigen::Index i = 0; i < mi; ++i) {
            d[static_cast<std::size_t>(i) + 1] = d_int(i);
        }

        // Compute smoothed positions: s = q - lambda * Q * d_int
        // @cite biagiotti2009 -- Sec. 4.4.5
        VecX const s_vec = q - lambda * Q * d_int;
        std::vector<Scalar> s(n_pts);
        for (Eigen::Index i = 0; i < ni; ++i) {
            s[static_cast<std::size_t>(i)] = s_vec(i);
        }

        // Compute cubic coefficients per span from smoothed positions and
        // second derivatives, using the standard cubic spline formula:
        //   q(t) = a + b*s + c*s^2 + d_coeff*s^3  where s = t - t_i
        // @cite biagiotti2009 -- Sec. 4.4, eq. (4.10)-(4.11)
        coeffs_.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            auto const si = s[i];
            auto const si1 = s[i + 1];
            auto const di = d[i];
            auto const di1 = d[i + 1];
            auto const hi = h[i];

            // Standard cubic spline coefficients from second derivatives (moments)
            // a = s_i
            // b = (s_{i+1} - s_i)/h_i - h_i*(2*d_i + d_{i+1})/6
            // c = d_i / 2
            // d_coeff = (d_{i+1} - d_i) / (6 * h_i)
            coeffs_[i][0] = si;
            coeffs_[i][1] = (si1 - si) / hi - hi * (Scalar{2} * di + di1) / Scalar{6};
            coeffs_[i][2] = di / Scalar{2};
            coeffs_[i][3] = (di1 - di) / (Scalar{6} * hi);
        }
    }

    /// @brief Find span index for time t using binary search.
    auto find_span(Scalar t) const -> std::size_t
    {
        auto const n_spans = coeffs_.size();
        auto it = std::upper_bound(times_.begin(), times_.end(), t);
        if (it == times_.begin()) {
            return 0;
        }
        auto idx = static_cast<std::size_t>(it - times_.begin()) - 1;
        if (idx >= n_spans) {
            idx = n_spans - 1;
        }
        return idx;
    }

    std::vector<Scalar> times_;
    std::vector<std::array<Scalar, 4>> coeffs_; ///< {a, b, c, d} per span
};

static_assert(trajectory_segment<smoothing_spline<double>, double, 1>);

} // namespace ctrlpp

#endif
