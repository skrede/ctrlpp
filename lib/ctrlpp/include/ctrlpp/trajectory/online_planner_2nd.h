#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_ONLINE_PLANNER_2ND_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_ONLINE_PLANNER_2ND_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2009, Sec. 4.6.2
// @cite lambrechts2005 -- Lambrechts, Boerlage & Steinbuch, "Trajectory Planning and Feedforward Design for Electromechanical Motion Systems", Control Engineering Practice 13(2), 2005

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/online_planner_diagnostics.h"

#include "ctrlpp/trajectory/detail/planner_2nd_evaluation.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include <cmath>
#include <algorithm>

namespace ctrlpp
{

// @cite biagiotti2009 -- Sec. 4.6.2
template <ctrlpp_floating_scalar Scalar>
class online_planner_2nd : private detail::planner_2nd_evaluation<Scalar>
{
    using base = detail::planner_2nd_evaluation<Scalar>;

    struct unchecked_t
    {
        explicit unchecked_t() = default;
    };

public:
    struct config
    {
        Scalar v_max;
        Scalar a_max;
        Scalar position_settle_tol = static_cast<Scalar>(1e-9);
        Scalar velocity_settle_tol = static_cast<Scalar>(1e-9);
    };

    // @cite biagiotti2009 -- Sec. 4.6.2
    static ctrlpp::expected<online_planner_2nd, trajectory_error> create(config const& cfg)
    {
        if(!std::isfinite(cfg.v_max) || cfg.v_max <= Scalar{0})
        {
            return ctrlpp::unexpected(trajectory_error::non_positive_velocity_limit);
        }
        if(!std::isfinite(cfg.a_max) || cfg.a_max <= Scalar{0})
        {
            return ctrlpp::unexpected(trajectory_error::non_positive_acceleration_limit);
        }
        return online_planner_2nd{unchecked_t{}, cfg};
    }

    // @cite biagiotti2009 -- Sec. 4.6.2
    ctrlpp::expected<void, trajectory_error> update(Scalar target)
    {
        if(!std::isfinite(target))
            return ctrlpp::unexpected(trajectory_error::non_finite_input);

        this->m_target = target;
        this->m_t_ref = this->m_t_last;

        this->m_q_ref = this->m_q;
        this->m_v_ref = this->m_v;

        this->compute_profile();
        return {};
    }

    // @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2009, Sec. 3.2.7, eq. (3.13a)-(3.13c), p.71
    trajectory_point<Scalar, 1> sample(Scalar t) const
    {
        auto const dt = std::max(t - this->m_t_ref, Scalar{0});

        Scalar q{};
        Scalar v{};
        Scalar a{};

        this->evaluate_profile(dt, q, v, a);

        this->m_q = q;
        this->m_v = v;
        this->m_t_last = t;

        this->m_settled = (std::abs(q - this->m_target) < this->m_position_settle_tol) && (std::abs(v) < this->m_velocity_settle_tol);

        return {
                .position = Vector<Scalar, 1>{q},
                .velocity = Vector<Scalar, 1>{v},
                .acceleration = Vector<Scalar, 1>{a},
        };
    }

    bool is_settled() const { return this->m_settled; }

    online_planner_diagnostics<Scalar> const& diagnostics() const { return this->m_diagnostics; }

    void reset(Scalar q0)
    {
        this->m_q = q0;
        this->m_v = Scalar{0};
        this->m_q_ref = q0;
        this->m_v_ref = Scalar{0};
        this->m_target = q0;
        this->m_t_ref = Scalar{0};
        this->m_t_last = Scalar{0};
        this->m_accel_duration = Scalar{0};
        this->m_cruise_duration = Scalar{0};
        this->m_decel_duration = Scalar{0};
        this->m_duration = Scalar{0};
        this->m_settled = true;
        this->record(online_planner_disposition::settled, online_planner_substitution_reason::none, Scalar{0}, Scalar{0}, q0);
    }

private:
    online_planner_2nd(unchecked_t, config const& cfg)
        : base{cfg}
    {
    }
};

}

#endif
