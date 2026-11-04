#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_SYNCHRONIZE_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_SYNCHRONIZE_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
// Automatic Machines and Robots", 2009, Sec. 5.3

#include "ctrlpp/trajectory/trajectory_types.h"

#include "ctrlpp/expected.h"

#include <span>
#include <concepts>
#include <algorithm>
#include <type_traits>

namespace ctrlpp
{

// @cite biagiotti2009 -- Sec. 5.3.1, eq. (5.15) -- duration-based rescaling preserves shape
template <typename T>
concept syncable_profile = requires(T& p, T const& cp, typename T::scalar_type dur) {
    { p.duration() } -> std::convertible_to<typename T::scalar_type>;
    { p.rescale_to(dur) } -> std::same_as<ctrlpp::expected<void, trajectory_error>>;
    { cp.can_rescale_to(dur) } -> std::same_as<ctrlpp::expected<void, trajectory_error>>;
};

// Every axis is checked before any is retimed, so a rejection leaves all of them unchanged. The slowest axis is
// retimed to a bit-exact copy of its own duration, which takes the no-op path.
// @cite biagiotti2009 -- Sec. 5.3, eq. (5.13)-(5.14) -- rescale_to concept for synchronization
template <syncable_profile... Profiles>
ctrlpp::expected<void, trajectory_error> synchronize(Profiles&... profiles)
{
    using Scalar = std::common_type_t<typename Profiles::scalar_type...>;
    auto const max_dur = std::max({static_cast<Scalar>(profiles.duration())...});

    auto status = ctrlpp::expected<void, trajectory_error>{};

    auto check = [&](auto& profile)
    {
        if(status.has_value())
            status = profile.can_rescale_to(max_dur);
    };
    (check(profiles), ...);
    if(!status.has_value())
        return status;

    auto commit = [&](auto& profile)
    {
        if(status.has_value())
            status = profile.rescale_to(max_dur);
    };
    (commit(profiles), ...);
    return status;
}

// @cite biagiotti2009 -- Sec. 5.3, eq. (5.14) -- time scaling for synchronization
template <syncable_profile Profile>
ctrlpp::expected<void, trajectory_error> synchronize(std::span<Profile> profiles)
{
    if(profiles.empty())
        return {};

    using Scalar = typename Profile::scalar_type;
    Scalar max_dur = profiles.front().duration();
    for(auto const& profile : profiles)
        max_dur = std::max(max_dur, static_cast<Scalar>(profile.duration()));

    for(auto const& profile : profiles)
    {
        auto const checked = profile.can_rescale_to(max_dur);
        if(!checked.has_value())
            return checked;
    }

    for(auto& profile : profiles)
    {
        auto const committed = profile.rescale_to(max_dur);
        if(!committed.has_value())
            return committed;
    }
    return {};
}

}

#endif
