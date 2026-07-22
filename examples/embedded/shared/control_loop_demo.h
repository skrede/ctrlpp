#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_CONTROL_LOOP_DEMO_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_CONTROL_LOOP_DEMO_H

#include "golden_reference.h"

#include "ctrlpp/control.h"

#include <Eigen/Dense>

#include <optional>

namespace ctrlpp
{

// Passive Scalar-parameterized control-loop kernel shared verbatim by every board
// leg: it owns the plant build and the on-device gain design, but no clock, no IO
// and no heap on step() -- the wrappers own the loop, cadence and telemetry.
template <class Scalar>
struct control_loop_demo
{
    Eigen::Matrix<Scalar, 2, 2> A;
    Eigen::Matrix<Scalar, 2, 1> B;
    Eigen::Matrix<Scalar, 1, 2> K;
    Eigen::Vector<Scalar, 2>    x;

    static std::optional<control_loop_demo> make()
    {
        const Scalar dt = static_cast<Scalar>(kDt);

        Eigen::Matrix<Scalar, 2, 2> A;
        A << Scalar{1}, dt, Scalar{0}, Scalar{1};
        Eigen::Matrix<Scalar, 2, 1> B;
        B << Scalar{0.5} * dt * dt, dt;
        Eigen::Matrix<Scalar, 2, 2> Q;
        Q << Scalar{10}, Scalar{0}, Scalar{0}, Scalar{1};
        Eigen::Matrix<Scalar, 1, 1> R;
        R << Scalar{0.1};

        const auto gain = ctrlpp::lqr_gain<Scalar, 2, 1>(A, B, Q, R);
        if(!gain.has_value())
            return std::nullopt;

        control_loop_demo demo;
        demo.A = A;
        demo.B = B;
        demo.K = *gain;
        demo.x << Scalar{1}, Scalar{0};
        return demo;
    }

    Scalar step()
    {
        const Scalar u = -(K * x)(0);
        x = A * x + B * u;
        return u;
    }
};

}

#endif
