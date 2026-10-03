#include "../../examples/embedded/shared/golden_verdict.h"
#include "../../examples/embedded/shared/control_loop_demo.h"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>

namespace {

// The Riccati count, written out here from the kernel's dimensions rather than
// by calling the header's helper, so a factor silently changed in the header
// stops matching.
template<class Scalar>
double recounted_gain_bound()
{
    const double m        = 2.0 * 2.0 + 1.0;
    const double assembly = 2.0 * 2.0 * 2.0 + 2.0 * 1.0 + 1.0 * 1.0;
    const double count    = assembly + 25.0 * m * m * m + 28.0 * m * m * m + 256.0 * m * m;
    const double norm     = std::abs(ctrlpp::kHostK0) + std::abs(ctrlpp::kHostK1);
    return count * static_cast<double>(std::numeric_limits<Scalar>::epsilon()) * norm;
}

template<class Scalar>
ctrlpp::control_loop_demo<Scalar> settled_run()
{
    auto demo = ctrlpp::control_loop_demo<Scalar>::make();
    REQUIRE(demo.has_value());
    for(int k = 0; k <= ctrlpp::kSteps; ++k)
        demo->step();
    return *demo;
}

}

TEST_CASE("The gain tolerance is the Riccati count recomputed from its factors", "[embedded][tolerance]")
{
    const double host = recounted_gain_bound<double>();
    REQUIRE(ctrlpp::golden_gain_tolerance<float>() == recounted_gain_bound<float>() + host);
    REQUIRE(ctrlpp::golden_gain_tolerance<double>() == host + host);
}

TEST_CASE("The float and double runs of the shared kernel pass beneath both bounds", "[embedded][tolerance]")
{
    ctrlpp::golden_bound bound          = ctrlpp::make_golden_bound();
    const ctrlpp::golden_verdict single = ctrlpp::judge_golden(settled_run<float>(), bound);
    const ctrlpp::golden_verdict twin   = ctrlpp::judge_golden(settled_run<double>(), bound);

    REQUIRE(single.pass);
    REQUIRE(twin.pass);
    WARN("float final-norm departure " << single.norm_departure << ", bound " << single.norm_tolerance << ", bound relative to the reference norm "
                                       << single.norm_tolerance / ctrlpp::kHostFinalNorm);
}

TEST_CASE("The float bound is tighter than the quantity it judges", "[embedded][tolerance]")
{
    ctrlpp::golden_bound bound          = ctrlpp::make_golden_bound();
    const ctrlpp::golden_verdict single = ctrlpp::judge_golden(settled_run<float>(), bound);

    REQUIRE(single.norm_tolerance < ctrlpp::kHostFinalNorm);
}

TEST_CASE("A board whose state settled to zero is refused", "[embedded][tolerance]")
{
    ctrlpp::golden_bound bound = ctrlpp::make_golden_bound();
    auto run                   = settled_run<float>();
    run.x.setZero();

    const ctrlpp::golden_verdict verdict = ctrlpp::judge_golden(run, bound);
    REQUIRE(verdict.gain_pass);
    REQUIRE_FALSE(verdict.norm_pass);
    REQUIRE_FALSE(verdict.pass);
}

TEST_CASE("A final state perturbed past its bound is refused and one within it accepted", "[embedded][tolerance]")
{
    ctrlpp::golden_bound bound = ctrlpp::make_golden_bound();
    const auto run             = settled_run<double>();
    const double reach         = ctrlpp::judge_golden(run, bound).norm_tolerance / ctrlpp::kHostFinalNorm;

    auto past = run;
    past.x *= 1.0 + 2.0 * reach;
    REQUIRE_FALSE(ctrlpp::judge_golden(past, bound).norm_pass);

    auto within = run;
    within.x *= 1.0 + 0.5 * reach;
    REQUIRE(ctrlpp::judge_golden(within, bound).norm_pass);
}

TEST_CASE("A gain perturbed past its bound is refused even when the state agrees", "[embedded][tolerance]")
{
    ctrlpp::golden_bound bound = ctrlpp::make_golden_bound();
    auto run                   = settled_run<double>();
    run.K(0, 0) += 2.0 * ctrlpp::golden_gain_tolerance<double>();

    const ctrlpp::golden_verdict verdict = ctrlpp::judge_golden(run, bound);
    REQUIRE_FALSE(verdict.gain_pass);
    REQUIRE_FALSE(verdict.pass);
}

TEST_CASE("The trajectory bound grows with the gain departure it is fed", "[embedded][tolerance]")
{
    ctrlpp::golden_bound bound = ctrlpp::make_golden_bound();
    const double eps           = static_cast<double>(std::numeric_limits<float>::epsilon());
    const double near          = bound.final_state_departure(eps, ctrlpp::kPlantRoundings, 1.0e-4);
    const double far           = bound.final_state_departure(eps, ctrlpp::kPlantRoundings, 2.0e-4);

    REQUIRE(far > near);
}

TEST_CASE("A departure exactly equal to its bound is refused", "[embedded][tolerance]")
{
    REQUIRE_FALSE(ctrlpp::within_bound(ctrlpp::golden_gain_tolerance<float>(), ctrlpp::golden_gain_tolerance<float>()));
    REQUIRE_FALSE(ctrlpp::within_bound(ctrlpp::kHostFinalNorm, ctrlpp::kHostFinalNorm));
}
