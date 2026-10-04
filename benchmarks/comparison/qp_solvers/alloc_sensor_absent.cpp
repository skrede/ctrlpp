#include "alloc_sensor.h"

namespace ctrlpp::bench::qp_solvers
{

auto alloc_sensor_present() -> bool
{
    return false;
}

void arm_alloc_sensor()
{
}

auto disarm_alloc_sensor() -> std::int64_t
{
    return -1;
}

}
