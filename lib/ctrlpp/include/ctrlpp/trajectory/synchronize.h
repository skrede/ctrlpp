#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_SYNCHRONIZE_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_SYNCHRONIZE_H

/// @brief Multi-axis trajectory synchronization.
///
/// Provides the synchronize() free function that rescales multiple trajectory
/// profiles to finish together at the duration of the slowest axis. Works with
/// any profile satisfying the syncable_profile concept, which requires a
/// duration() query, a fallible rescale_to(), and the matching fallible
/// can_rescale_to() query that makes the all-or-nothing guarantee below possible.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
/// Automatic Machines and Robots", 2009, Sec. 5.3

#include "ctrlpp/expected.h"

#include "ctrlpp/trajectory/trajectory_types.h"

#include <span>
#include <concepts>
#include <algorithm>
#include <type_traits>

namespace ctrlpp
{

/// @brief Concept for profiles that support duration query and time rescaling.
///
/// A syncable profile exposes a scalar_type alias, a duration() method returning
/// a value convertible to scalar_type, a rescale_to() method that returns the
/// fallible result of retiming the profile, and a const can_rescale_to() query
/// with the same return that answers the same question without mutating. The
/// rescaling return is constrained exactly, so a profile whose rescaling cannot
/// report failure fails at the concept rather than at the call site.
///
/// @cite biagiotti2009 -- Sec. 5.3
/// @cite biagiotti2009 -- Sec. 5.3.1, eq. (5.15) -- duration-based rescaling preserves shape
template <typename T>
concept syncable_profile = requires(T& p, T const& cp, typename T::scalar_type dur) {
    { p.duration() } -> std::convertible_to<typename T::scalar_type>;
    { p.rescale_to(dur) } -> std::same_as<ctrlpp::expected<void, trajectory_error>>;
    { cp.can_rescale_to(dur) } -> std::same_as<ctrlpp::expected<void, trajectory_error>>;
};

/// @brief Synchronize multiple axes to finish simultaneously.
///
/// Finds the maximum duration across all profiles and retimes each to it. The
/// work runs in two passes: every axis is asked whether it can reach the target
/// duration, and only once all of them have answered yes is any axis retimed. A
/// rejection on any axis returns that axis's error with nothing mutated. Both
/// passes go through the same solve inside each profile, so the checking pass
/// cannot pass an axis the committing pass then fails, and neither pass builds an
/// owning copy of anything.
///
/// The slowest axis never trips the shortening rejection: the target is a
/// bit-exact copy of that axis's own reported duration, so it compares equal and
/// takes the success no-op path. There is no float-equality fragility in that,
/// because the value compared is the axis's own stored duration rather than a
/// recomputed one.
///
/// @cite biagiotti2009 -- Sec. 5.3, eq. (5.13)-(5.14) -- rescale_to concept for synchronization
template <syncable_profile... Profiles>
[[nodiscard]] auto synchronize(Profiles&... profiles) -> ctrlpp::expected<void, trajectory_error>
{
    using Scalar = std::common_type_t<typename Profiles::scalar_type...>;
    auto const max_dur = std::max({static_cast<Scalar>(profiles.duration())...});

    auto status = ctrlpp::expected<void, trajectory_error>{};

    auto check = [&](auto& profile) {
        if (status.has_value()) {
            status = profile.can_rescale_to(max_dur);
        }
    };
    (check(profiles), ...);
    if (!status.has_value()) {
        return status;
    }

    auto commit = [&](auto& profile) {
        if (status.has_value()) {
            status = profile.rescale_to(max_dur);
        }
    };
    (commit(profiles), ...);
    return status;
}

/// @brief Synchronize a contiguous run of homogeneous profiles.
///
/// Runtime-sized variant. The view carries no ownership, so the same call serves
/// a vector, a plain array, or a statically allocated buffer and nothing on this
/// path allocates. Construct one from any contiguous container at the call site,
/// for example `ctrlpp::synchronize(std::span{axes})`.
///
/// Two passes, all-or-nothing, exactly as above.
///
/// @cite biagiotti2009 -- Sec. 5.3, eq. (5.14) -- time scaling for synchronization
template <syncable_profile Profile>
[[nodiscard]] auto synchronize(std::span<Profile> profiles) -> ctrlpp::expected<void, trajectory_error>
{
    if (profiles.empty()) {
        return {};
    }

    using Scalar = typename Profile::scalar_type;
    Scalar max_dur = profiles.front().duration();
    for (auto const& profile : profiles) {
        max_dur = std::max(max_dur, static_cast<Scalar>(profile.duration()));
    }

    for (auto const& profile : profiles) {
        auto const checked = profile.can_rescale_to(max_dur);
        if (!checked.has_value()) {
            return checked;
        }
    }

    for (auto& profile : profiles) {
        auto const committed = profile.rescale_to(max_dur);
        if (!committed.has_value()) {
            return committed;
        }
    }
    return {};
}

}

#endif
