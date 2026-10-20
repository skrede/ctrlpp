#include "ctrlpp/estimation/particle_filter.h"

#include <random>

struct canary_dynamics
{
    auto operator()(const ctrlpp::Vector<double, 1>& x, const ctrlpp::Vector<double, 1>&) const -> ctrlpp::Vector<double, 1> { return x; }
};

struct canary_measurement
{
    auto operator()(const ctrlpp::Vector<double, 1>& x) const -> ctrlpp::Vector<double, 1> { return x; }
};

using canary_filter = ctrlpp::particle_filter<double, 1, 1, 1, 4, canary_dynamics, canary_measurement>;

int main()
{
    const ctrlpp::pf_config<double, 1, 1, 1> cfg{};
#ifdef CTRLPP_CANARY_PF_DIRECT_CONSTRUCTION
    canary_filter filter{{}, canary_dynamics{}, canary_measurement{}, cfg, {}, ctrlpp::systematic_resampling{}, std::mt19937_64{}, ctrlpp::detail::gaussian_likelihood<double, 1>{}};
    return filter.particles()[0].allFinite() ? 0 : 1;
#else
    const auto filter = canary_filter::create(canary_dynamics{}, canary_measurement{}, cfg);
    return filter.has_value() ? 0 : 1;
#endif
}
