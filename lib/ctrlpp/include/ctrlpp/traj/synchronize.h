#ifndef HPP_GUARD_CTRLPP_TRAJ_SYNCHRONIZE_H
#define HPP_GUARD_CTRLPP_TRAJ_SYNCHRONIZE_H

/// @brief Multi-axis trajectory synchronization.
///
/// Provides the synchronize() free function that rescales multiple trajectory
/// profiles to finish simultaneously at the duration of the slowest axis.
/// Works with any profile satisfying the syncable_profile concept (must have
/// duration() and rescale_to() methods).
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
/// Automatic Machines and Robots", 2009, Sec. 5.3

#include <algorithm>
#include <concepts>
#include <type_traits>
#include <vector>

namespace ctrlpp
{

/// @brief Concept for profiles that support duration query and time rescaling.
///
/// A syncable profile must expose a scalar_type alias, a duration() method
/// returning a value convertible to scalar_type, and a rescale_to() method
/// that accepts a new target duration.
///
/// @cite biagiotti2009 -- Sec. 5.3
/// @cite biagiotti2009 -- Sec. 5.3.1, eq. (5.15) -- duration-based rescaling preserves shape
template <typename T>
concept syncable_profile = requires(T& p, typename T::scalar_type dur) {
    { p.duration() } -> std::convertible_to<typename T::scalar_type>;
    { p.rescale_to(dur) };
};

/// @brief Synchronize multiple axes to finish simultaneously.
///
/// Finds the maximum duration across all profiles and rescales each to that
/// duration. Each profile's rescale_to() is responsible for maintaining
/// constraint satisfaction (v_max, a_max, j_max).
///
/// @cite biagiotti2009 -- Sec. 5.3, eq. (5.13)-(5.14) -- rescale_to concept for synchronization
template <syncable_profile... Profiles>
void synchronize(Profiles&... profiles)
{
    using Scalar = std::common_type_t<typename Profiles::scalar_type...>;
    auto const max_dur = std::max({static_cast<Scalar>(profiles.duration())...});
    (profiles.rescale_to(max_dur), ...);
}

/// @brief Synchronize a vector of homogeneous profiles to finish simultaneously.
///
/// Runtime-sized variant for collections of identical profile types.
/// Uses rescale_to() which preserves kinematic limit satisfaction.
///
/// @cite biagiotti2009 -- Sec. 5.3, eq. (5.14) -- time scaling for synchronization
template <syncable_profile Profile>
void synchronize(std::vector<Profile>& profiles)
{
    if (profiles.empty()) {
        return;
    }
    using Scalar = typename Profile::scalar_type;
    Scalar max_dur{};
    for (auto const& p : profiles) {
        max_dur = std::max(max_dur, p.duration());
    }
    for (auto& p : profiles) {
        p.rescale_to(max_dur);
    }
}

} // namespace ctrlpp

#endif
