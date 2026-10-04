#include "driver_options.h"

#include <array>
#include <ranges>
#include <cstdio>
#include <utility>
#include <charconv>
#include <algorithm>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

auto split(std::string_view text) -> std::vector<std::string>
{
    std::vector<std::string> parts;
    for(const auto part : std::views::split(text, ','))
        if(!part.empty())
            parts.emplace_back(part.begin(), part.end());
    return parts;
}

auto integer(std::string_view text, std::int32_t& out) -> bool
{
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), out);
    return error == std::errc{} && end == text.data() + text.size();
}

using list_member = std::vector<std::string> driver_options::*;
using text_member = std::string driver_options::*;

constexpr std::array<std::pair<std::string_view, list_member>, 6> list_options{{{"--cell", &driver_options::cells},
                                                                                {"--variant", &driver_options::variants},
                                                                                {"--regime", &driver_options::regimes},
                                                                                {"--tier", &driver_options::tiers},
                                                                                {"--arm", &driver_options::arms},
                                                                                {"--refused-families", &driver_options::refused_families}}};

constexpr std::array<std::pair<std::string_view, text_member>, 5> text_options{{{"--out", &driver_options::out_dir},
                                                                                {"--run-id", &driver_options::run_id},
                                                                                {"--describe", &driver_options::describe_path},
                                                                                {"--binding-gate", &driver_options::gate_path},
                                                                                {"--on-probe-failure", &driver_options::on_probe_failure}}};

auto apply(driver_options& options, std::string_view key, std::string_view value) -> bool
{
    for(const auto& [name, member] : list_options)
        if(key == name)
            return !(options.*member = split(value)).empty();
    for(const auto& [name, member] : text_options)
        if(key == name)
            return !(options.*member = value).empty();
    if(key == "--rep")
        return integer(value, options.rep) && options.rep >= 0 && options.rep < repetitions;
    if(key == "--steps")
        return integer(value, options.steps) && options.steps > 0 && options.steps <= steps_per_sequence;
    return key == "--cohort" && value == name_of(options.group);
}

}

auto parse_options(cohort group, int argc, char** argv) -> std::optional<driver_options>
{
    driver_options options{group, ".", "", 0, steps_per_sequence, {}, {}, {}, {}, {}, {}, "", "", "abort"};
    for(int i = 1; i + 1 < argc; i += 2)
        if(!apply(options, argv[i], argv[i + 1]))
        {
            std::fprintf(stderr, "qp_solvers_%s: bad argument %s %s\n", name_of(group).data(), argv[i], argv[i + 1]);
            return std::nullopt;
        }
    if(options.on_probe_failure != "abort" && options.on_probe_failure != "continue")
    {
        std::fprintf(stderr, "qp_solvers_%s: --on-probe-failure is abort or continue\n", name_of(group).data());
        return std::nullopt;
    }
    if(argc % 2 == 0)
    {
        std::fprintf(stderr, "qp_solvers_%s: every option takes one value\n", name_of(group).data());
        return std::nullopt;
    }
    return options;
}

auto admits(const std::vector<std::string>& filter, std::string_view name) -> bool
{
    return filter.empty() || std::ranges::find(filter, name) != filter.end();
}

}
