#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_ALLOC_SENSOR_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_ALLOC_SENSOR_H

#include <cstdint>

namespace ctrlpp {

struct canary_verdict
{
    bool live;
    std::uint32_t observed;
    std::uint32_t expected;
    std::uint32_t eigen_trips;
};

// Counting is off until armed, so setup allocation is never charged to the
// window being measured. Arming also forbids Eigen allocation, which is what
// makes the Eigen sentinel trip.
void alloc_sensor_arm() noexcept;

void alloc_sensor_disarm() noexcept;

void alloc_sensor_reset() noexcept;

std::uint32_t alloc_sensor_observed() noexcept;

std::uint32_t alloc_sensor_eigen_trips() noexcept;

// Makes one deliberate allocation through each path the sensor claims to see
// inside an armed window. A zero the sensor reports later is evidence only if
// this came back live.
canary_verdict run_alloc_canary() noexcept;

}

#endif
