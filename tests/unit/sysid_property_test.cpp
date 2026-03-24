#include "ctrlpp/sysid/batch_arx.h"
#include "ctrlpp/sysid/rls.h"

#include <catch2/catch_test_macros.hpp>
#include <rapidcheck.h>
#include <rapidcheck/catch.h>

#include <cmath>
#include <cstddef>

using namespace ctrlpp;

namespace
{

auto bounded_double(double lo, double hi) -> rc::Gen<double>
{
    return rc::gen::map(rc::gen::inRange(0, 1000000), [lo, hi](int x) { return lo + (hi - lo) * (static_cast<double>(x) / 1000000.0); });
}

} // namespace

TEST_CASE("sysid property tests", "[sysid][property]")
{
    SECTION("rls identified model order matches configured order")
    {
        rc::prop("rls model order consistency", [](void)
                 {
            constexpr std::size_t NP = 4;
            rls<double, NP> identifier;

            // Feed random data points
            for(int i = 0; i < 50; ++i)
            {
                auto y_val = *bounded_double(-5.0, 5.0);
                Vector<double, NP> phi;
                for(std::size_t j = 0; j < NP; ++j)
                    phi(static_cast<Eigen::Index>(j)) = *bounded_double(-5.0, 5.0);
                identifier.update(y_val, phi);
            }

            // Check that parameters vector has correct dimension
            auto params = identifier.parameters();
            RC_ASSERT(params.size() == static_cast<Eigen::Index>(NP));
            for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NP); ++i)
                RC_ASSERT(std::isfinite(params(i))); });
    }

    SECTION("batch_arx model order consistency")
    {
        rc::prop("batch_arx order", [](void)
                 {
            constexpr std::size_t NA = 2;
            constexpr std::size_t NB = 2;
            constexpr std::size_t N_SAMPLES = 50;

            // Generate random signals
            Eigen::Matrix<double, 1, Eigen::Dynamic> Y(1, static_cast<Eigen::Index>(N_SAMPLES));
            Eigen::Matrix<double, 1, Eigen::Dynamic> U(1, static_cast<Eigen::Index>(N_SAMPLES));

            for(std::size_t i = 0; i < N_SAMPLES; ++i)
            {
                Y(0, static_cast<Eigen::Index>(i)) = *bounded_double(-5.0, 5.0);
                U(0, static_cast<Eigen::Index>(i)) = *bounded_double(-5.0, 5.0);
            }

            auto result = batch_arx<NA, NB>(Y, U);

            // State-space system A should be NA x NA, B should be NA x 1
            RC_ASSERT(result.system.A.rows() == static_cast<Eigen::Index>(NA));
            RC_ASSERT(result.system.A.cols() == static_cast<Eigen::Index>(NA));
            RC_ASSERT(result.system.B.rows() == static_cast<Eigen::Index>(NA));
            RC_ASSERT(result.system.B.cols() == 1);

            // All system matrices should be finite
            RC_ASSERT(std::isfinite(result.system.A.norm()));
            RC_ASSERT(std::isfinite(result.system.B.norm())); });
    }
}
