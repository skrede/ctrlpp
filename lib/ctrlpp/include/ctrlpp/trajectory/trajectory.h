#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_TRAJECTORY_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
// Machines and Robots", 2009, Sec. 5.2.1, eq. (2.16)-(2.17), p.34

#include "ctrlpp/trajectory/cycloidal_path.h"
#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/trajectory_segment.h"

#include "ctrlpp/expected.h"

#include <cmath>
#include <cstddef>
#include <utility>
#include <algorithm>

namespace ctrlpp
{

template <typename Law, typename Scalar, std::size_t ND>
class trajectory
{
public:
    trajectory(Law law, Vector<Scalar, ND> q0, Vector<Scalar, ND> displacement, Scalar duration)
        : m_law{std::move(law)}
        , m_q0{q0}
        , m_displacement{displacement}
        , m_duration{duration}
    {
    }

    trajectory_point<Scalar, ND> evaluate(Scalar t) const
    {
        auto const tau = std::clamp(t, Scalar{0}, m_duration) / m_duration;
        auto const ml = m_law(tau);
        auto const inv_T = Scalar{1} / m_duration;
        return {
                .position = m_q0 + m_displacement * ml.q,
                .velocity = m_displacement * (inv_T * ml.dq),
                .acceleration = m_displacement * (inv_T * inv_T * ml.ddq),
        };
    }

    Scalar duration() const { return m_duration; }

private:
    Law m_law;
    Vector<Scalar, ND> m_q0;
    Vector<Scalar, ND> m_displacement;
    Scalar m_duration;
};

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
//   Machines and Robots", 2009, Ch. 2-5 -- complete trajectory generation framework
template <typename Law, typename Scalar, int Rows>
ctrlpp::expected<trajectory<Law, Scalar, static_cast<std::size_t>(Rows)>, trajectory_error> make_trajectory(
        Law law, Eigen::Matrix<Scalar, Rows, 1> const& q0, Eigen::Matrix<Scalar, Rows, 1> const& q1, Scalar duration)
{
    constexpr auto ND = static_cast<std::size_t>(Rows);
    if(!std::isfinite(duration) || !q0.allFinite() || !q1.allFinite())
        return ctrlpp::unexpected(trajectory_error::non_finite_input);
    if(duration <= Scalar{0})
        return ctrlpp::unexpected(trajectory_error::non_positive_duration);
    Vector<Scalar, ND> const displacement = (q1 - q0).eval();
    return trajectory<Law, Scalar, ND>{std::move(law), q0, displacement, duration};
}

static_assert(trajectory_segment<trajectory<path_point<double> (*)(double), double, 1>, double, 1>);

}

#endif
