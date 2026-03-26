#include <ctrlpp/traj/cubic_path.h>
#include <ctrlpp/traj/cubic_trajectory.h>
#include <ctrlpp/traj/cycloidal_path.h>
#include <ctrlpp/traj/piecewise_path.h>
#include <ctrlpp/traj/piecewise_trajectory.h>
#include <ctrlpp/traj/quintic_path.h>
#include <ctrlpp/traj/quintic_trajectory.h>
#include <ctrlpp/traj/trajectory.h>
#include <ctrlpp/traj/trajectory_segment.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace
{

auto zero1() -> Eigen::Matrix<double, 1, 1> { return Eigen::Matrix<double, 1, 1>::Zero(); }
auto val1(double v) -> Eigen::Matrix<double, 1, 1> { return Eigen::Matrix<double, 1, 1>::Constant(v); }

} // namespace

// ---- piecewise_trajectory tests ----

TEST_CASE("piecewise_trajectory with two cubic trajectories", "[traj][piecewise]")
{
    auto seg1 = ctrlpp::make_cubic_trajectory(val1(0.0), val1(5.0), zero1(), zero1(), 1.0);
    auto seg2 = ctrlpp::make_cubic_trajectory(val1(5.0), val1(10.0), zero1(), zero1(), 1.0);
    auto pw = ctrlpp::piecewise_trajectory<double, 1, decltype(seg1), decltype(seg2)>{seg1, seg2};

    SECTION("duration is sum of segment durations")
    {
        REQUIRE_THAT(pw.duration(), WithinAbs(2.0, 1e-12));
    }

    SECTION("evaluate(0) matches first segment start")
    {
        auto pt = pw.evaluate(0.0);
        REQUIRE_THAT(pt.position[0], WithinAbs(0.0, 1e-12));
    }

    SECTION("evaluate(0.999) matches first segment at t=0.999")
    {
        auto expected = seg1.evaluate(0.999);
        auto actual = pw.evaluate(0.999);
        REQUIRE_THAT(actual.position[0], WithinAbs(expected.position[0], 1e-10));
    }

    SECTION("evaluate(1.0) matches second segment at t=0 (boundary)")
    {
        auto expected = seg2.evaluate(0.0);
        auto actual = pw.evaluate(1.0);
        REQUIRE_THAT(actual.position[0], WithinAbs(expected.position[0], 1e-10));
        REQUIRE_THAT(actual.position[0], WithinAbs(5.0, 1e-10));
    }

    SECTION("evaluate(2.0) matches second segment endpoint")
    {
        auto pt = pw.evaluate(2.0);
        REQUIRE_THAT(pt.position[0], WithinAbs(10.0, 1e-10));
    }

    SECTION("evaluate(0.5) is between 0 and 5")
    {
        auto pt = pw.evaluate(0.5);
        REQUIRE(pt.position[0] > 0.0);
        REQUIRE(pt.position[0] < 5.0);
    }
}

TEST_CASE("piecewise_trajectory with cubic + quintic (heterogeneous)", "[traj][piecewise]")
{
    auto cubic = ctrlpp::make_cubic_trajectory(val1(0.0), val1(10.0), zero1(), zero1(), 2.0);
    auto quintic = ctrlpp::make_quintic_trajectory(val1(10.0), val1(20.0), zero1(), zero1(), zero1(), zero1(), 3.0);

    auto pw = ctrlpp::piecewise_trajectory<double, 1, decltype(cubic), decltype(quintic)>{cubic, quintic};

    SECTION("duration is sum")
    {
        REQUIRE_THAT(pw.duration(), WithinAbs(5.0, 1e-12));
    }

    SECTION("evaluate(0) is start")
    {
        REQUIRE_THAT(pw.evaluate(0.0).position[0], WithinAbs(0.0, 1e-10));
    }

    SECTION("evaluate(2) is boundary")
    {
        REQUIRE_THAT(pw.evaluate(2.0).position[0], WithinAbs(10.0, 1e-10));
    }

    SECTION("evaluate(5) is end")
    {
        REQUIRE_THAT(pw.evaluate(5.0).position[0], WithinAbs(20.0, 1e-10));
    }
}

TEST_CASE("piecewise_trajectory with three segments: cubic + quintic + trajectory cycloidal", "[traj][piecewise]")
{
    auto cubic = ctrlpp::make_cubic_trajectory(val1(0.0), val1(5.0), zero1(), zero1(), 1.0);
    auto quintic = ctrlpp::make_quintic_trajectory(val1(5.0), val1(15.0), zero1(), zero1(), zero1(), zero1(), 2.0);
    auto cycloidal = ctrlpp::make_trajectory(ctrlpp::cycloidal_path<double>, val1(15.0), val1(20.0), 1.5);

    auto pw = ctrlpp::piecewise_trajectory<double, 1, decltype(cubic), decltype(quintic), decltype(cycloidal)>{
        cubic, quintic, cycloidal};

    SECTION("duration is sum of three")
    {
        REQUIRE_THAT(pw.duration(), WithinAbs(4.5, 1e-12));
    }

    SECTION("evaluates at boundaries")
    {
        REQUIRE_THAT(pw.evaluate(0.0).position[0], WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(pw.evaluate(1.0).position[0], WithinAbs(5.0, 1e-10));
        REQUIRE_THAT(pw.evaluate(3.0).position[0], WithinAbs(15.0, 1e-10));
        REQUIRE_THAT(pw.evaluate(4.5).position[0], WithinAbs(20.0, 1e-10));
    }
}

TEST_CASE("piecewise_trajectory clamping", "[traj][piecewise]")
{
    auto seg1 = ctrlpp::make_cubic_trajectory(val1(0.0), val1(5.0), zero1(), zero1(), 1.0);
    auto seg2 = ctrlpp::make_cubic_trajectory(val1(5.0), val1(10.0), zero1(), zero1(), 1.0);
    auto pw = ctrlpp::piecewise_trajectory<double, 1, decltype(seg1), decltype(seg2)>{seg1, seg2};

    SECTION("evaluate before 0 returns first segment start")
    {
        auto at_neg = pw.evaluate(-1.0);
        auto at_zero = pw.evaluate(0.0);
        REQUIRE_THAT(at_neg.position[0], WithinAbs(at_zero.position[0], 1e-12));
    }

    SECTION("evaluate past duration returns last segment endpoint")
    {
        auto at_end = pw.evaluate(pw.duration());
        auto at_past = pw.evaluate(100.0);
        REQUIRE_THAT(at_past.position[0], WithinAbs(at_end.position[0], 1e-12));
    }
}

TEST_CASE("piecewise_trajectory satisfies trajectory_segment concept", "[traj][piecewise]")
{
    using cubic_t = ctrlpp::cubic_trajectory<double, 1>;
    using quintic_t = ctrlpp::quintic_trajectory<double, 1>;
    using pw_t = ctrlpp::piecewise_trajectory<double, 1, cubic_t, quintic_t>;
    static_assert(ctrlpp::trajectory_segment<pw_t, double, 1>);
    SUCCEED("piecewise_trajectory satisfies trajectory_segment concept");
}

TEST_CASE("piecewise_trajectory continuity at joints", "[traj][piecewise]")
{
    auto seg1 = ctrlpp::make_cubic_trajectory(val1(0.0), val1(5.0), zero1(), zero1(), 1.0);
    auto seg2 = ctrlpp::make_cubic_trajectory(val1(5.0), val1(10.0), zero1(), zero1(), 1.0);
    auto pw = ctrlpp::piecewise_trajectory<double, 1, decltype(seg1), decltype(seg2)>{seg1, seg2};

    constexpr double eps = 1e-12;
    auto before = pw.evaluate(1.0 - eps);
    auto at = pw.evaluate(1.0);

    REQUIRE_THAT(before.position[0], WithinAbs(at.position[0], 1e-6));
    REQUIRE_THAT(before.position[0], WithinAbs(5.0, 1e-6));
    REQUIRE_THAT(at.position[0], WithinAbs(5.0, 1e-10));
}

// ---- piecewise_path tests ----

namespace
{

/// @brief Simple wrapper turning a path function into a path_segment for piecewise_path testing.
template <typename Scalar, ctrlpp::path_point<Scalar> (*Fn)(Scalar)>
struct path_segment_wrapper
{
    Scalar dur;
    auto evaluate(Scalar tau) const -> ctrlpp::path_point<Scalar> { return Fn(tau / dur); }
    auto duration() const -> Scalar { return dur; }
};

} // namespace

TEST_CASE("piecewise_path composes two cubic paths", "[traj][piecewise_path]")
{
    using cubic_seg = path_segment_wrapper<double, ctrlpp::cubic_path<double>>;
    auto seg1 = cubic_seg{.dur = 1.0};
    auto seg2 = cubic_seg{.dur = 1.0};

    auto pw = ctrlpp::piecewise_path<double, cubic_seg, cubic_seg>{seg1, seg2};

    SECTION("duration is sum")
    {
        REQUIRE_THAT(pw.duration(), WithinAbs(2.0, 1e-12));
    }

    SECTION("evaluate(0) gives q=0")
    {
        auto pt = pw.evaluate(0.0);
        REQUIRE_THAT(pt.q, WithinAbs(0.0, 1e-12));
    }

    SECTION("evaluate(1.0) gives q=1 (end of first segment)")
    {
        auto pt = pw.evaluate(1.0);
        REQUIRE_THAT(pt.q, WithinAbs(0.0, 1e-10));
    }

    SECTION("evaluate(2.0) gives q=1 (end of second segment)")
    {
        auto pt = pw.evaluate(2.0);
        REQUIRE_THAT(pt.q, WithinAbs(1.0, 1e-10));
    }
}

TEST_CASE("piecewise_path composes cubic + quintic (heterogeneous)", "[traj][piecewise_path]")
{
    using cubic_seg = path_segment_wrapper<double, ctrlpp::cubic_path<double>>;
    using quintic_seg = path_segment_wrapper<double, ctrlpp::quintic_path<double>>;
    auto seg1 = cubic_seg{.dur = 2.0};
    auto seg2 = quintic_seg{.dur = 3.0};

    auto pw = ctrlpp::piecewise_path<double, cubic_seg, quintic_seg>{seg1, seg2};

    SECTION("duration is sum")
    {
        REQUIRE_THAT(pw.duration(), WithinAbs(5.0, 1e-12));
    }

    SECTION("midpoint of first segment has q in (0, 1)")
    {
        auto pt = pw.evaluate(1.0);
        REQUIRE(pt.q > 0.0);
        REQUIRE(pt.q < 1.0);
    }
}
