#include "ctrlpp/control/care.h"
#include "ctrlpp/control/dare.h"

#include "ctrlpp/detail/riccati_solution.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Core>

#include <cmath>
#include <random>
#include <numbers>

namespace
{

using matrix3 = Eigen::Matrix<double, 3, 3>;
using vector3 = Eigen::Matrix<double, 3, 1>;

constexpr int rotations = 64;

struct rotated_pose
{
    matrix3 A;
    vector3 B;
    matrix3 Q;
    matrix3 P;
};

// Raw engine output mapped by explicit arithmetic, because the standard distributions differ between standard libraries.
auto unit_interval(std::mt19937_64 &engine) -> double
{
    return static_cast<double>(engine() >> 11) * 0x1p-53;
}

// Rodrigues' rotation formula, T = I + sin(t) K + (1 - cos(t)) K^2 with K the cross-product matrix of a unit axis.
auto rodrigues_rotation(std::mt19937_64 &engine) -> matrix3
{
    vector3 axis;
    for(int i = 0; i < 3; ++i)
        axis(i) = 2.0 * unit_interval(engine) - 1.0;
    axis.normalize();
    const double angle = std::numbers::pi * (2.0 * unit_interval(engine) - 1.0);
    matrix3 K;
    K << 0.0, -axis(2), axis(1), axis(2), 0.0, -axis(0), -axis(1), axis(0), 0.0;
    return matrix3::Identity() + std::sin(angle) * K + (1.0 - std::cos(angle)) * K * K;
}

// The base pose drives and weights the first mode alone, so its solution is diag(p, 0, 0) with p the scalar Riccati root.
auto rotated(const matrix3 &T, const vector3 &modes, double p) -> rotated_pose
{
    const vector3 e1 = vector3::UnitX();
    const matrix3 Q0 = e1 * e1.transpose();
    const matrix3 TQT = T * Q0 * T.transpose();
    return {T * modes.asDiagonal() * T.transpose(), T * e1, (TQT + TQT.transpose()) / 2, T * (p * Q0) * T.transpose()};
}

auto relative_distance(const matrix3 &P, const matrix3 &reference) -> double
{
    return (P - reference).norm() / reference.norm();
}

template<typename Solve, typename Refusal>
auto accepted_rotations(const vector3 &modes, double p, Solve solve, Refusal is_expected_refusal) -> int
{
    std::mt19937_64 engine;
    int accepted = 0;
    for(int k = 0; k < rotations; ++k)
    {
        const rotated_pose pose = rotated(rodrigues_rotation(engine), modes, p);
        const auto result       = solve(pose);
        CAPTURE(k);
        if(result)
        {
            ++accepted;
            CHECK(relative_distance(result->P, pose.P) <= ctrlpp::detail::half_significand_margin<double>());
        }
        else
        {
            CHECK(is_expected_refusal(result.error()));
        }
    }
    return accepted;
}

const Eigen::Matrix<double, 1, 1> unit_weight = Eigen::Matrix<double, 1, 1>::Identity();

}

TEST_CASE("A rotated rank-one continuous regulator is accepted wherever its answer is semi-definite to rounding", "[riccati][definiteness]")
{
    const vector3 modes{1.0, -1.0, -2.0};
    const double p = 1.0 + std::sqrt(2.0);
    const int accepted = accepted_rotations(
        modes, p, [](const rotated_pose &pose) { return ctrlpp::care<double, 3, 1>(pose.A, pose.B, pose.Q, unit_weight); },
        [](ctrlpp::care_error error) { return error == ctrlpp::care_error::non_psd_solution; });
    CAPTURE(accepted);
    // Each rotation's exact solution is semi-definite, so a refusal here would decline a correct answer.
    CHECK(accepted == rotations);
}

TEST_CASE("A rotated rank-one discrete regulator is accepted wherever its answer is semi-definite to rounding", "[riccati][definiteness]")
{
    const double a = 1.5;
    const vector3 modes{a, 0.5, -0.25};
    const double p = (a * a + std::sqrt(a * a * a * a + 4.0)) / 2;
    const int accepted = accepted_rotations(
        modes, p, [](const rotated_pose &pose) { return ctrlpp::dare<double, 3, 1>(pose.A, pose.B, pose.Q, unit_weight); },
        [](ctrlpp::dare_error error) { return error == ctrlpp::dare_error::non_psd_solution || error == ctrlpp::dare_error::arithmetic_limit; });
    CAPTURE(accepted);
    // Each rotation's exact solution is semi-definite, so a refusal here would decline a correct answer.
    CHECK(accepted == rotations);
}
