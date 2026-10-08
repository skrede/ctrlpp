#include "build_id.h"
#include "telemetry.h"
#include "sbrk_ceiling.h"

#include "golden_reference.h"

#include <cstdio>
#include <cstdint>
#include <cinttypes>

namespace ctrlpp {

namespace {

constexpr std::uint64_t kTenthsPerCycle = 10;

static_assert(column_count(claim_class::spread) == 11, "emit_spread_row prints eleven values in the table's order");

void print_view(std::string_view text)
{
    std::printf("%.*s", static_cast<int>(text.size()), text.data());
}

void print_bits(double value)
{
    char text[kBitsTextLength<double>];
    std::printf(",%s", render_bits(value, text).data());
}

}

void open_line(stream_tag tag)
{
    print_view(tag_text(tag));
}

void emit_identity()
{
    char build_id[kBuildIdTextLength];
    render_build_id(build_id);
    open_line(stream_tag::meta);
    std::printf("build_id=%s\n", build_id);
    open_line(stream_tag::meta);
    std::printf("kind=%s\n", kAllocSensorPresent ? "evidence" : "timing");
    open_line(stream_tag::meta);
    std::printf("board=NUCLEO-H753ZI loop=bare-superloop scalar=double\n");
}

void emit_allocation_canary(const canary_verdict &verdict)
{
    open_line(stream_tag::canary);
    std::printf("allocation %s observed=%" PRIu32 " expected=%" PRIu32 " eigen_trips=%" PRIu32 "\n", verdict.live ? "PASS" : "FAIL", verdict.observed, verdict.expected,
                verdict.eigen_trips);
}

void emit_counter_canary(counter_state state)
{
    open_line(stream_tag::canary);
    if(state == counter_state::live)
        std::printf("cycle_counter PASS state=live\n");
    else
        std::printf("cycle_counter UNAVAILABLE -- %s; no timing figure is printed\n", describe(state));
}

void emit_header(claim_class claim)
{
    open_line(tag_of(claim));
    const char *separator = "";
    for(const claim_field &field : kClaimFields)
    {
        if(field.claim != claim)
            continue;
        std::printf("%s", separator);
        print_view(field.name);
        separator = ",";
    }
    std::printf("\n");
}

void emit_exact_row(std::int32_t k, const std::array<double, kExactValues> &values)
{
    open_line(stream_tag::exact);
    std::printf("%" PRId32, k);
    for(const double value : values)
        print_bits(value);
    std::printf("\n");
}

void emit_spread_row(const timing_region &region, const cycle_distribution &d)
{
    open_line(stream_tag::spread);
    std::printf("%s,%s,%" PRIu32 ",%s,%" PRIu32, region.name, region.regime, d.n, region.basis, region.block);
    if(d.n == 0)
    {
        std::printf(",,,,,,%" PRIu32 "\n", d.wraps);
        return;
    }
    const std::uint64_t tenths = (std::uint64_t{d.max} * kTenthsPerCycle + region.block / 2U) / region.block;
    std::printf(",%" PRIu32 ",%" PRIu32 ",%" PRIu32 ",%" PRIu32, d.min, d.median, d.max, d.max_at);
    std::printf(",%" PRIu32 ".%" PRIu32 ",%" PRIu32 "\n", static_cast<std::uint32_t>(tenths / kTenthsPerCycle), static_cast<std::uint32_t>(tenths % kTenthsPerCycle), d.wraps);
}

void emit_setup_heap()
{
    open_line(stream_tag::verdict);
    std::printf("heap setup_high_water=%" PRIu32 " reserve=%" PRIu32 " bytes\n", static_cast<std::uint32_t>(heap_high_water_bytes()), static_cast<std::uint32_t>(heap_reserve_bytes()));
}

void emit_gain(double k0, double k1)
{
    open_line(stream_tag::verdict);
    std::printf("gain device=%.17g,%.17g host=%.17g,%.17g\n", k0, k1, kHostK0, kHostK1);
}

}
