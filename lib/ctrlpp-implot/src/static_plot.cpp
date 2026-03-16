#include "ctrlpp/implot/static_plot.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

namespace ctrlpp::implot {

namespace {

constexpr int margin_l = 70;
constexpr int margin_r = 20;
constexpr int margin_t = 30;
constexpr int margin_b = 50;

constexpr const char* palette[] = {
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"
};
constexpr int palette_size = 6;

struct DataRange {
    float t_min;
    float t_max;
    float y_min;
    float y_max;
};

struct ChronoPoint {
    float t;
    float v;
};

auto read_chronological(const RingBuffer& buf) -> std::vector<ChronoPoint>
{
    std::vector<ChronoPoint> pts;
    pts.reserve(static_cast<std::size_t>(buf.size));
    for (int i = 0; i < buf.size; ++i) {
        int idx = (buf.offset - buf.size + i + buf.max_size) % buf.max_size;
        pts.push_back({buf.data[idx].x, buf.data[idx].y});
    }
    return pts;
}

auto compute_range(const SignalRecorder& recorder,
                   std::initializer_list<const char*> signals) -> std::optional<DataRange>
{
    float t_min = std::numeric_limits<float>::max();
    float t_max = std::numeric_limits<float>::lowest();
    float y_min = std::numeric_limits<float>::max();
    float y_max = std::numeric_limits<float>::lowest();
    bool has_data = false;

    for (const auto* name : signals) {
        const auto* buf = recorder.channel(name);
        if (!buf || buf->size == 0) continue;
        has_data = true;
        for (int i = 0; i < buf->size; ++i) {
            int idx = (buf->offset - buf->size + i + buf->max_size) % buf->max_size;
            float t = buf->data[idx].x;
            float v = buf->data[idx].y;
            t_min = std::min(t_min, t);
            t_max = std::max(t_max, t);
            y_min = std::min(y_min, v);
            y_max = std::max(y_max, v);
        }
    }

    if (!has_data) return std::nullopt;

    // Pad y range
    if (y_min == y_max) {
        y_min -= 1.0f;
        y_max += 1.0f;
    } else {
        float pad = (y_max - y_min) * 0.05f;
        y_min -= pad;
        y_max += pad;
    }

    return DataRange{t_min, t_max, y_min, y_max};
}

void write_tick_labels(std::ofstream& out, int plot_w, int plot_h,
                       const DataRange& range)
{
    // Y-axis ticks
    for (int i = 0; i < 5; ++i) {
        float frac = static_cast<float>(i) / 4.0f;
        float val = range.y_min + frac * (range.y_max - range.y_min);
        float py = static_cast<float>(plot_h) * (1.0f - frac);
        out << "    <line x1=\"-4\" y1=\"" << std::fixed << std::setprecision(2) << py
            << "\" x2=\"0\" y2=\"" << py
            << "\" stroke=\"#333\" stroke-width=\"1\"/>\n";
        out << "    <text x=\"-8\" y=\"" << py + 4.0f
            << "\" font-size=\"11\" fill=\"#333\" text-anchor=\"end\">"
            << std::setprecision(2) << val << "</text>\n";
    }

    // X-axis ticks
    for (int i = 0; i < 5; ++i) {
        float frac = static_cast<float>(i) / 4.0f;
        float val = range.t_min + frac * (range.t_max - range.t_min);
        float px = static_cast<float>(plot_w) * frac;
        out << "    <line x1=\"" << std::fixed << std::setprecision(2) << px
            << "\" y1=\"" << plot_h
            << "\" x2=\"" << px
            << "\" y2=\"" << plot_h + 4
            << "\" stroke=\"#333\" stroke-width=\"1\"/>\n";
        out << "    <text x=\"" << px
            << "\" y=\"" << plot_h + 18
            << "\" font-size=\"11\" fill=\"#333\" text-anchor=\"middle\">"
            << std::setprecision(2) << val << "</text>\n";
    }
}

}

void write_svg(const char* path,
               const SignalRecorder& recorder,
               std::initializer_list<const char*> signals,
               const SvgConfig& config)
{
    auto range = compute_range(recorder, signals);
    if (!range) return;

    int plot_w = config.width - margin_l - margin_r;
    int plot_h = config.height - margin_t - margin_b;

    std::ofstream out(path);
    if (!out) return;

    out << std::fixed << std::setprecision(2);

    // SVG header
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 "
        << config.width << " " << config.height << "\">\n"
        << "  <rect width=\"" << config.width << "\" height=\""
        << config.height << "\" fill=\"white\"/>\n"
        << "  <g transform=\"translate(" << margin_l << "," << margin_t << ")\">\n";

    // Axes
    out << "    <line x1=\"0\" y1=\"0\" x2=\"0\" y2=\"" << plot_h
        << "\" stroke=\"#333\" stroke-width=\"1\"/>\n";
    out << "    <line x1=\"0\" y1=\"" << plot_h << "\" x2=\"" << plot_w
        << "\" y2=\"" << plot_h << "\" stroke=\"#333\" stroke-width=\"1\"/>\n";

    // Tick marks and labels
    write_tick_labels(out, plot_w, plot_h, *range);

    // Title
    out << "    <text x=\"" << plot_w / 2
        << "\" y=\"-10\" font-size=\"14\" font-weight=\"bold\" fill=\"#333\" text-anchor=\"middle\">"
        << config.title << "</text>\n";

    // X-axis label
    out << "    <text x=\"" << plot_w / 2
        << "\" y=\"" << plot_h + 40
        << "\" font-size=\"12\" fill=\"#333\" text-anchor=\"middle\">"
        << config.x_label << "</text>\n";

    // Y-axis label (rotated)
    out << "    <text x=\"-55\" y=\"" << plot_h / 2
        << "\" font-size=\"12\" fill=\"#333\" text-anchor=\"middle\""
        << " transform=\"rotate(-90,-55," << plot_h / 2 << ")\">"
        << config.y_label << "</text>\n";

    // Signal polylines
    float t_span = range->t_max - range->t_min;
    float y_span = range->y_max - range->y_min;
    if (t_span == 0.0f) t_span = 1.0f;
    if (y_span == 0.0f) y_span = 1.0f;

    int color_idx = 0;
    int legend_y = 10;

    for (const auto* name : signals) {
        const auto* buf = recorder.channel(name);
        if (!buf || buf->size == 0) {
            ++color_idx;
            continue;
        }

        const char* color = palette[color_idx % palette_size];
        auto pts = read_chronological(*buf);

        out << "    <polyline fill=\"none\" stroke=\"" << color
            << "\" stroke-width=\"1.5\" points=\"";

        for (std::size_t i = 0; i < pts.size(); ++i) {
            float px = static_cast<float>(plot_w) * (pts[i].t - range->t_min) / t_span;
            float py = static_cast<float>(plot_h) * (1.0f - (pts[i].v - range->y_min) / y_span);
            if (i > 0) out << " ";
            out << px << "," << py;
        }

        out << "\"/>\n";

        // Legend entry
        int lx = plot_w - 120;
        out << "    <line x1=\"" << lx << "\" y1=\"" << legend_y
            << "\" x2=\"" << lx + 20 << "\" y2=\"" << legend_y
            << "\" stroke=\"" << color << "\" stroke-width=\"2\"/>\n";
        out << "    <text x=\"" << lx + 25 << "\" y=\"" << legend_y + 4
            << "\" font-size=\"11\" fill=\"#333\">" << name << "</text>\n";

        legend_y += 18;
        ++color_idx;
    }

    out << "  </g>\n</svg>\n";
}

}
