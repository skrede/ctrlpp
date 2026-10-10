#include "online_planner_3rd_case.h"

#include <cstddef>
#include <cstdint>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    return ctrlpp::fuzz::run_online_planner_3rd_case<float>(data, size);
}
