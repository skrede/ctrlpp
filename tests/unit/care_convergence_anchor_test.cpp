#include "ctrlpp/control/care.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <random>
#include <cstddef>
#include <cstdint>
#include <numbers>

namespace
{

constexpr std::uint64_t sweep_seed = 0xc0a57ab1e5ULL;
constexpr int decade_count = 18;
constexpr int cases_per_decade = 64;
constexpr int cases_per_scale_direction = 32;

template <typename Scalar>
using matrix2 = Eigen::Matrix<Scalar, 2, 2>;

template <typename Scalar>
auto closed_loop_is_stable(const matrix2<Scalar>& A,
                           const matrix2<Scalar>& B,
                           const matrix2<Scalar>& R,
                           const matrix2<Scalar>& P) -> bool
{
    const matrix2<Scalar> closed_loop =
        (A - B * R.inverse() * B.transpose() * P).eval();
    Eigen::EigenSolver<matrix2<Scalar>> eigensystem(closed_loop, false);
    for(int index = 0; index < 2; ++index)
    {
        if(!(eigensystem.eigenvalues()(index).real() < Scalar{0}))
            return false;
    }
    return true;
}

// Weight scale up to which the solver is required to carry the posing. It is a
// conservative floor, not the measured boundary: weights within a factor
// 1/epsilon of unit dynamics stay separable from them in double arithmetic, so
// nothing in that band may be declined. A fine log-grid scan of both families
// below puts the true accept/decline boundary at 1e16.20 (comfortable) and
// 1e16.40 (degenerate), leaving roughly half a decade of margin above this
// floor -- enough that a different `uniform_real_distribution` sequence cannot
// move a draw across it.
const double representable_weight_scale =
    1.0 / std::numeric_limits<double>::epsilon();

auto is_enumerated_care_error(ctrlpp::care_error error) -> bool
{
    switch(error)
    {
    case ctrlpp::care_error::non_lhp_stabilizable:
    case ctrlpp::care_error::non_finite_input:
    case ctrlpp::care_error::singular_r:
    case ctrlpp::care_error::singular_u11:
    case ctrlpp::care_error::non_psd_solution:
    case ctrlpp::care_error::schur_failed:
    case ctrlpp::care_error::sign_function_stagnated:
        return true;
    }
    return false;
}

auto rotation(double angle) -> matrix2<double>
{
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    matrix2<double> result;
    result << cosine, -sine, sine, cosine;
    return result;
}

}

TEST_CASE("CARE sign iteration declines an unresolved regenerated input",
          "[care][convergence][error]")
{
    matrix2<double> A;
    A << 0x0p+0, 0x0p+0,
        -0xf.dfdfdfdfdfdf8p-7, -0xf.dfdfdfdfdfdf8p-7;

    Eigen::Matrix<double, 2, 1> B;
    B << -0xf.dfdfdfdfdfdf8p-7,
         -0xf.dfdfdfdfdfdf8p-7;

    matrix2<double> raw_weight;
    raw_weight << -0xf.dfdfdfdfffdf8p-7, -0xf.dfdfdfdfdfdf8p-7,
                  -0xf.dfdfdfdfdf98p-7,  -0xf.dfdfdfdfdfdf8p-7;
    const matrix2<double> Q =
        (raw_weight.transpose() * raw_weight).eval();

    Eigen::Matrix<double, 1, 1> R;
    R << 0x8.3126e978d4fep-11;

    const auto result = ctrlpp::care<double, 2, 1>(A, B, Q, R);

    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::care_error::sign_function_stagnated);

    const auto reference = ctrlpp::care<long double, 2, 1>(
        A.cast<long double>(),
        B.cast<long double>(),
        Q.cast<long double>(),
        R.cast<long double>());
    CHECK_FALSE(reference.has_value());
}

TEST_CASE("CARE sign iteration accepts only stable near-axis solutions",
          "[care][convergence][sweep]")
{
    std::mt19937_64 generator(sweep_seed);
    std::uniform_real_distribution<double> mantissa(1.0, 10.0);
    std::uniform_real_distribution<double> alpha_distribution(0.2, 0.95);
    std::uniform_real_distribution<double> angle_distribution(
        -std::numbers::pi, std::numbers::pi);

    std::size_t near_axis_drawn = 0;
    std::size_t near_axis_accepted = 0;
    std::size_t near_axis_accepted_unstable = 0;

    // Acceptance is asserted per draw, against the draw's own weight scale,
    // rather than per decade against a count. A per-decade count is a property
    // of the generator's draw sequence, not of the solver: libstdc++, libc++
    // and MSVC give `std::uniform_real_distribution` different sequences from
    // the same seed, so an equality on the count fails on two of the three
    // while the solver behaves identically.
    std::size_t representable_drawn = 0;
    std::size_t representable_accepted = 0;
    std::size_t extreme_drawn = 0;
    std::size_t extreme_accepted = 0;

    for(int decade = 1; decade <= decade_count; ++decade)
    {
        const double decade_scale = std::pow(10.0, -decade);
        for(int index = 0; index < cases_per_decade; ++index)
        {
            const double delta = decade_scale * mantissa(generator);
            const double alpha = alpha_distribution(generator);
            const double beta = std::sqrt(1.0 - alpha * alpha);
            const double angle = angle_distribution(generator);
            const matrix2<double> transform = rotation(angle);

            matrix2<double> diagonal_a;
            diagonal_a << alpha * delta, 0.0, 0.0, -0.8;
            matrix2<double> diagonal_b;
            diagonal_b << beta * delta, 0.0, 0.0, 0.6;
            const matrix2<double> A =
                (transform * diagonal_a * transform.transpose()).eval();
            const matrix2<double> B =
                (transform * diagonal_b).eval();
            const matrix2<double> Q = matrix2<double>::Identity();
            const matrix2<double> R = matrix2<double>::Identity();

            const auto result = ctrlpp::care<double, 2, 2>(A, B, Q, R);
            ++near_axis_drawn;
            if(result.has_value())
            {
                ++near_axis_accepted;
                near_axis_accepted_unstable +=
                    !closed_loop_is_stable(A, B, R, result->P);
                CAPTURE(decade, index, delta, alpha, angle);
                CHECK(std::isnan(result->subspace_separation));
                CHECK(result->reorder_complete);
            }
        }

        for(int direction : {-1, 1})
        {
            for(int index = 0;
                index < cases_per_scale_direction;
                ++index)
            {
                const double exponent =
                    static_cast<double>(direction * decade);
                const double scale =
                    std::pow(10.0, exponent) * mantissa(generator);
                const matrix2<double> transform =
                    rotation(angle_distribution(generator));

                matrix2<double> base_a;
                base_a << -1.1, 0.3, -0.2, -1.4;
                matrix2<double> base_b;
                base_b << 0.8, 0.1, -0.15, 0.65;
                matrix2<double> base_q;
                base_q << 1.0, 0.2, 0.2, 1.7;
                matrix2<double> base_r;
                base_r << 1.3, 0.1, 0.1, 0.9;
                const matrix2<double> comfortable_a =
                    (transform * base_a * transform.transpose()).eval();
                const matrix2<double> comfortable_b =
                    (transform * base_b * transform.transpose()).eval();
                const matrix2<double> comfortable_q =
                    (scale * transform * base_q * transform.transpose()).eval();
                const matrix2<double> comfortable_r =
                    (scale * transform * base_r * transform.transpose()).eval();
                const bool representable = scale <= representable_weight_scale;
                if(representable)
                    representable_drawn += 2;
                else
                    extreme_drawn += 2;

                const auto comfortable = ctrlpp::care<double, 2, 2>(
                    comfortable_a,
                    comfortable_b,
                    comfortable_q,
                    comfortable_r);
                CAPTURE(decade, direction, index, scale, representable);
                if(comfortable.has_value())
                {
                    if(representable)
                        ++representable_accepted;
                    else
                        ++extreme_accepted;
                    CHECK(closed_loop_is_stable(
                        comfortable_a,
                        comfortable_b,
                        comfortable_r,
                        comfortable->P));
                }
                else
                {
                    // Beyond the guaranteed band a decline is permitted, but it
                    // must still reach the caller as an enumerated cause.
                    CHECK(is_enumerated_care_error(comfortable.error()));
                }

                const matrix2<double> degenerate_a =
                    -matrix2<double>::Identity();
                const matrix2<double> degenerate_b =
                    matrix2<double>::Identity();
                const matrix2<double> degenerate_q =
                    scale * matrix2<double>::Identity();
                const matrix2<double> degenerate_r =
                    scale * matrix2<double>::Identity();
                const auto degenerate = ctrlpp::care<double, 2, 2>(
                    degenerate_a,
                    degenerate_b,
                    degenerate_q,
                    degenerate_r);
                if(degenerate.has_value())
                {
                    if(representable)
                        ++representable_accepted;
                    else
                        ++extreme_accepted;
                    CHECK(closed_loop_is_stable(
                        degenerate_a,
                        degenerate_b,
                        degenerate_r,
                        degenerate->P));
                }
                else
                {
                    CHECK(is_enumerated_care_error(degenerate.error()));
                }
            }
        }
    }

    CAPTURE(sweep_seed,
            near_axis_drawn,
            near_axis_accepted,
            near_axis_accepted_unstable);
    CHECK(near_axis_drawn
          == static_cast<std::size_t>(decade_count * cases_per_decade));
    CHECK(near_axis_accepted > 0);
    CHECK(near_axis_accepted_unstable == 0);

    // Inside the guaranteed band the solver must not decline anything: this is
    // the over-rejection guard, and it is an equality because every draw there
    // is required to succeed. Outside it, acceptance is a property of the
    // problem's conditioning and is deliberately not asserted -- the inline
    // checks above still require every accepted answer to be stabilizing and
    // every decline to be enumerated.
    CAPTURE(sweep_seed,
            representable_weight_scale,
            representable_drawn,
            representable_accepted,
            extreme_drawn,
            extreme_accepted);
    CHECK(representable_drawn > 0);
    CHECK(extreme_drawn > 0);
    CHECK(representable_accepted == representable_drawn);
}

TEST_CASE("CARE sign iteration preserves its diagnostic contract",
          "[care][convergence][diagnostic]")
{
    const matrix2<float> A = -matrix2<float>::Identity();
    const matrix2<float> B = matrix2<float>::Identity();
    const matrix2<float> Q = matrix2<float>::Identity();
    const matrix2<float> R = matrix2<float>::Identity();

    const auto result = ctrlpp::care<float, 2, 2>(A, B, Q, R);

    REQUIRE(result.has_value());
    CHECK(closed_loop_is_stable(A, B, R, result->P));
    CHECK(std::isnan(result->subspace_separation));
    CHECK(result->reorder_complete);
}
