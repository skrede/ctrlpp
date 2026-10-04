#include "frequency_log.h"

#include <format>

#if defined(__linux__)
    #include <sched.h>
#endif

namespace ctrlpp::bench::qp_solvers
{

frequency_log::frequency_log(const std::string& path)
    : out_(path, std::ios::app)
{
}

void frequency_log::sample(std::string_view label)
{
#if defined(__linux__)
    const int cpu = sched_getcpu();
    std::ifstream source(std::format("/sys/devices/system/cpu/cpu{}/cpufreq/scaling_cur_freq", cpu));
    long long khz = 0;
    if(cpu >= 0 && source >> khz)
        out_ << label << ' ' << cpu << ' ' << khz << '\n';
#else
    static_cast<void>(label);
#endif
}

}
