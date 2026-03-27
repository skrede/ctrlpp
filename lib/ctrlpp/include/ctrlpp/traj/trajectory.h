#ifndef HPP_GUARD_CTRLPP_TRAJ_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJ_TRAJECTORY_H

/// @brief Adapter lifting normalized paths into physical trajectory segments.
///
/// Given a callable law(tau) -> path_point producing normalized values in [0,1],
/// trajectory applies kinematic scaling to produce physical position, velocity,
/// and acceleration for a point-to-point motion with displacement h over duration T.
///
/// Derivative scaling: vel = h/T * dq, acc = h/T^2 * ddq
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009, Sec. 5.2.1, eq. (2.16)-(2.17), p.34

#include "ctrlpp/traj/cycloidal_path.h"
#include "ctrlpp/traj/trajectory_segment.h"
#include "ctrlpp/traj/trajectory_types.h"

#include <algorithm>
#include <cstddef>
#include <utility>

namespace ctrlpp
{

/// @brief Adapter that lifts a normalized path into a trajectory_segment.
///
/// The law must be a callable with signature (Scalar tau) -> path_point<Scalar>.
template <typename Law, typename Scalar, std::size_t ND>
class trajectory
{
public:
    trajectory(Law law, Vector<Scalar, ND> q0, Vector<Scalar, ND> displacement, Scalar duration)
        : law_{std::move(law)}, q0_{q0}, h_{displacement}, dur_{duration}
    {
    }

    auto evaluate(Scalar t) const -> trajectory_point<Scalar, ND>
    {
        auto const tau = std::clamp(t, Scalar{0}, dur_) / dur_;
        auto const ml = law_(tau);
        auto const inv_T = Scalar{1} / dur_;
        return {
            .position = q0_ + h_ * ml.q,
            .velocity = h_ * (inv_T * ml.dq),
            .acceleration = h_ * (inv_T * inv_T * ml.ddq),
        };
    }

    auto duration() const -> Scalar { return dur_; }

private:
    Law law_;
    Vector<Scalar, ND> q0_;
    Vector<Scalar, ND> h_;
    Scalar dur_;
};

/// @brief Create a trajectory from a path, start/end positions, and duration.
///
/// Internally stores displacement h = q1 - q0.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
///   Machines and Robots", 2009, Ch. 2-5 -- complete trajectory generation framework
template <typename Law, typename Scalar, int Rows>
auto make_trajectory(
    Law law,
    Eigen::Matrix<Scalar, Rows, 1> const& q0,
    Eigen::Matrix<Scalar, Rows, 1> const& q1,
    Scalar duration) -> trajectory<Law, Scalar, static_cast<std::size_t>(Rows)>
{
    constexpr auto ND = static_cast<std::size_t>(Rows);
    Vector<Scalar, ND> const displacement = (q1 - q0).eval();
    return {std::move(law), q0, displacement, duration};
}

/// Verify trajectory satisfies trajectory_segment concept.
static_assert(trajectory_segment<
    trajectory<path_point<double> (*)(double), double, 1>, double, 1>);

} // namespace ctrlpp

#endif
