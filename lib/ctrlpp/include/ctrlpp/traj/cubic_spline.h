#ifndef HPP_GUARD_CTRLPP_TRAJ_CUBIC_SPLINE_H
#define HPP_GUARD_CTRLPP_TRAJ_CUBIC_SPLINE_H

/// @brief Cubic spline interpolation with natural, clamped, and periodic BCs.
///
/// Constructs a C2-continuous cubic spline through waypoints using the
/// velocity-based formulation. Three boundary condition types are supported:
/// natural (zero endpoint acceleration), clamped (assigned endpoint velocities),
/// and periodic (cyclic with matching start/end derivatives).
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
/// Automatic Machines and Robots", 2009, Sec. 4.4, eq. (4.10)-(4.11)

#include "ctrlpp/traj/detail/tridiagonal.h"
#include "ctrlpp/traj/trajectory_segment.h"
#include "ctrlpp/traj/trajectory_types.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <vector>

namespace ctrlpp
{

/// @brief Boundary condition type for cubic spline interpolation.
enum class boundary_condition
{
    natural,   ///< Zero acceleration at endpoints (M_0 = M_n = 0)
    clamped,   ///< Assigned endpoint velocities (v_0 and v_n specified)
    periodic   ///< Cyclic: v_0 = v_n and a_0 = a_n (requires q_0 = q_n)
};

/// @brief Cubic spline interpolation through waypoints.
///
/// Evaluates position, velocity, and acceleration at arbitrary time t using
/// piecewise cubic polynomials with C2 continuity at interior knots.
///
/// @cite biagiotti2009 -- Sec. 4.4, eq. (4.10)-(4.11)
template <typename Scalar>
class cubic_spline
{
  public:
    struct config
    {
        std::vector<Scalar> times;       ///< Knot times t_0 ... t_n (n+1 entries)
        std::vector<Scalar> positions;   ///< Waypoint positions q_0 ... q_n (n+1 entries)
        boundary_condition bc{boundary_condition::natural};
        Scalar v0{};                     ///< Endpoint velocity for clamped BC
        Scalar vn{};                     ///< Endpoint velocity for clamped BC
    };

    /// @brief Construct cubic spline from waypoints and boundary conditions.
    ///
    /// Sets up and solves the tridiagonal system for spline velocities,
    /// then computes polynomial coefficients for each span.
    ///
    /// @cite biagiotti2009 -- Sec. 4.4, eq. (4.10)-(4.11)
    explicit cubic_spline(config const& cfg)
        : times_{cfg.times}
    {
        auto const n_pts = times_.size();
        assert(n_pts >= 2);
        assert(cfg.positions.size() == n_pts);

        auto const n = n_pts - 1; // number of spans

        // Compute span durations h_i and slopes delta_i
        std::vector<Scalar> h(n);
        std::vector<Scalar> delta(n);
        for (std::size_t i = 0; i < n; ++i) {
            h[i] = times_[i + 1] - times_[i];
            assert(h[i] > Scalar{0});
            delta[i] = (cfg.positions[i + 1] - cfg.positions[i]) / h[i];
        }

        // Solve for velocities at each waypoint
        std::vector<Scalar> vel(n_pts);

        if (cfg.bc == boundary_condition::natural) {
            solve_natural(h, delta, cfg.positions, vel, n_pts, n);
        } else if (cfg.bc == boundary_condition::clamped) {
            solve_clamped(h, delta, cfg.positions, vel, n_pts, n, cfg.v0, cfg.vn);
        } else {
            assert(std::abs(cfg.positions.front() - cfg.positions.back()) < Scalar{1e-10});
            solve_periodic(h, delta, cfg.positions, vel, n_pts, n);
        }

        // Compute polynomial coefficients for each span
        // q(s) = a + b*s + c*s^2 + d*s^3  where s = t - t_i
        // @cite biagiotti2009 -- Sec. 4.4, eq. (4.10)-(4.11)
        coeffs_.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            auto const qi = cfg.positions[i];
            auto const qi1 = cfg.positions[i + 1];
            auto const vi = vel[i];
            auto const vi1 = vel[i + 1];
            auto const hi = h[i];

            coeffs_[i][0] = qi;                                                    // a
            coeffs_[i][1] = vi;                                                    // b
            coeffs_[i][2] = (Scalar{3} * (qi1 - qi) / hi - Scalar{2} * vi - vi1) / hi;  // c
            coeffs_[i][3] = (Scalar{2} * (qi - qi1) / hi + vi + vi1) / (hi * hi);       // d
        }
    }

    /// @brief Evaluate spline at time t, clamped to [t_0, t_n].
    ///
    /// Uses binary search to find the active span, then Horner evaluation
    /// for position, velocity, and acceleration.
    ///
    /// @cite biagiotti2009 -- Sec. 4.4
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
    /// @brief Find span index for time t using binary search.
    ///
    /// Clamps to valid span range [0, n_spans-1].
    auto find_span(Scalar t) const -> std::size_t
    {
        auto const n_spans = coeffs_.size();
        // upper_bound gives first element > t
        auto it = std::upper_bound(times_.begin(), times_.end(), t);
        if (it == times_.begin()) {
            return 0;
        }
        auto idx = static_cast<std::size_t>(it - times_.begin()) - 1;
        // Clamp to valid range (Pitfall 6)
        if (idx >= n_spans) {
            idx = n_spans - 1;
        }
        return idx;
    }

    /// @brief Solve natural BC: full (n+1)x(n+1) tridiagonal for all velocities.
    ///
    /// Row 0: 2*h_0 * v_0 + h_0 * v_1 = 3*(q_1 - q_0)  (from M_0=0)
    /// Row n: h_{n-1} * v_{n-1} + 2*h_{n-1} * v_n = 3*(q_n - q_{n-1})  (from M_n=0)
    /// Interior: h_i * v_{i-1} + 2*(h_{i-1}+h_i) * v_i + h_{i-1} * v_{i+1} = RHS_i
    ///
    /// @cite biagiotti2009 -- Sec. 4.4, eq. (4.10)-(4.11)
    static void solve_natural(std::vector<Scalar> const& h,
                              std::vector<Scalar> const& delta,
                              std::vector<Scalar> const& /*pos*/,
                              std::vector<Scalar>& vel,
                              std::size_t n_pts,
                              std::size_t n)
    {
        // Full system of size n+1
        std::vector<Scalar> a(n_pts, Scalar{0});
        std::vector<Scalar> b(n_pts, Scalar{0});
        std::vector<Scalar> c(n_pts, Scalar{0});
        std::vector<Scalar> d(n_pts, Scalar{0});

        // Row 0: natural BC (M_0 = 0 => 2*h_0*v_0 + h_0*v_1 = 3*(q_1 - q_0))
        b[0] = Scalar{2} * h[0];
        c[0] = h[0];
        d[0] = Scalar{3} * delta[0] * h[0];

        // Interior rows
        for (std::size_t i = 1; i < n; ++i) {
            a[i] = h[i];
            b[i] = Scalar{2} * (h[i - 1] + h[i]);
            c[i] = h[i - 1];
            d[i] = Scalar{3} * (h[i] * delta[i - 1] + h[i - 1] * delta[i]);
        }

        // Row n: natural BC (M_n = 0 => h_{n-1}*v_{n-1} + 2*h_{n-1}*v_n = 3*(q_n - q_{n-1}))
        a[n] = h[n - 1];
        b[n] = Scalar{2} * h[n - 1];
        d[n] = Scalar{3} * delta[n - 1] * h[n - 1];

        detail::thomas_solve(a, b, c, d);

        for (std::size_t i = 0; i < n_pts; ++i) {
            vel[i] = d[i];
        }
    }

    /// @brief Solve clamped BC: (n-1)x(n-1) system for interior velocities.
    ///
    /// v_0 and v_n are known from config. Interior system moves known terms to RHS.
    ///
    /// @cite biagiotti2009 -- Sec. 4.4.1
    static void solve_clamped(std::vector<Scalar> const& h,
                              std::vector<Scalar> const& delta,
                              std::vector<Scalar> const& /*pos*/,
                              std::vector<Scalar>& vel,
                              std::size_t /*n_pts*/,
                              std::size_t n,
                              Scalar v0,
                              Scalar vn)
    {
        vel[0] = v0;
        vel[n] = vn;

        if (n <= 1) {
            // Only 2 waypoints -- no interior unknowns
            return;
        }

        auto const n_int = n - 1; // number of interior unknowns (v_1 .. v_{n-1})
        std::vector<Scalar> a(n_int, Scalar{0});
        std::vector<Scalar> b(n_int, Scalar{0});
        std::vector<Scalar> c(n_int, Scalar{0});
        std::vector<Scalar> d(n_int, Scalar{0});

        for (std::size_t j = 0; j < n_int; ++j) {
            auto const i = j + 1; // original index
            a[j] = (j > 0) ? h[i] : Scalar{0};
            b[j] = Scalar{2} * (h[i - 1] + h[i]);
            c[j] = (j < n_int - 1) ? h[i - 1] : Scalar{0};
            d[j] = Scalar{3} * (h[i] * delta[i - 1] + h[i - 1] * delta[i]);

            // Move known endpoint velocities to RHS
            if (i == 1) {
                d[j] -= h[i] * v0;
            }
            if (i == n - 1) {
                d[j] -= h[i - 1] * vn;
            }
        }

        detail::thomas_solve(a, b, c, d);

        for (std::size_t j = 0; j < n_int; ++j) {
            vel[j + 1] = d[j];
        }
    }

    /// @brief Solve periodic BC: n x n cyclic tridiagonal for v_0 .. v_{n-1}.
    ///
    /// Requires q_0 == q_n. Velocity v_n = v_0 (periodic wrap-around).
    /// Uses cyclic Thomas solver (Sherman-Morrison).
    ///
    /// @cite biagiotti2009 -- Sec. 4.4.2
    static void solve_periodic(std::vector<Scalar> const& h,
                               std::vector<Scalar> const& delta,
                               std::vector<Scalar> const& /*pos*/,
                               std::vector<Scalar>& vel,
                               std::size_t /*n_pts*/,
                               std::size_t n)
    {
        // System of size n for v_0 .. v_{n-1} (v_n = v_0)
        std::vector<Scalar> a(n, Scalar{0});
        std::vector<Scalar> b(n, Scalar{0});
        std::vector<Scalar> c(n, Scalar{0});
        std::vector<Scalar> d(n, Scalar{0});

        // Row 0: wraps around -- v_{n-1} appears as sub-diagonal corner
        b[0] = Scalar{2} * (h[n - 1] + h[0]);
        c[0] = h[n - 1];
        d[0] = Scalar{3} * (h[0] * delta[n - 1] + h[n - 1] * delta[0]);

        // Interior rows 1..n-2
        for (std::size_t i = 1; i < n - 1; ++i) {
            a[i] = h[i];
            b[i] = Scalar{2} * (h[i - 1] + h[i]);
            c[i] = h[i - 1];
            d[i] = Scalar{3} * (h[i] * delta[i - 1] + h[i - 1] * delta[i]);
        }

        // Row n-1: wraps around -- v_0 appears as super-diagonal corner
        a[n - 1] = h[n - 1];
        b[n - 1] = Scalar{2} * (h[n - 2] + h[n - 1]);
        d[n - 1] = Scalar{3} * (h[n - 1] * delta[n - 2] + h[n - 2] * delta[n - 1]);

        // Corner elements: alpha = A[0][n-1] = h[0], beta = A[n-1][0] = h[n-2]
        auto const alpha = h[0];
        auto const beta = h[n - 2];

        detail::cyclic_thomas_solve(a, b, c, d, alpha, beta);

        for (std::size_t i = 0; i < n; ++i) {
            vel[i] = d[i];
        }
        vel[n] = vel[0]; // periodic wrap
    }

    std::vector<Scalar> times_;
    std::vector<std::array<Scalar, 4>> coeffs_; ///< {a, b, c, d} per span
};

static_assert(trajectory_segment<cubic_spline<double>, double, 1>);

} // namespace ctrlpp

#endif
