#include "ctrlpp/implot/scrolling_plot.h"

#include "implot.h"

namespace ctrlpp::implot {

static_assert(sizeof(ctrlpp::implot::Point2f) == 2 * sizeof(float));

void scrolling_plot(const PlotConfig& config, const SignalRecorder& recorder,
                    std::initializer_list<const char*> signal_names,
                    float current_time)
{
    if (ImPlot::BeginPlot(config.title, ImVec2(-1, config.height))) {
        ImPlot::SetupAxes(config.x_label, config.y_label);
        ImPlot::SetupAxisLimits(ImAxis_X1,
                                static_cast<double>(current_time - config.history),
                                static_cast<double>(current_time),
                                ImGuiCond_Always);

        for (const auto* name : signal_names) {
            const auto* buf = recorder.channel(name);
            if (buf && buf->size > 0) {
                ImPlot::PlotLine(name, &buf->data[0].first, &buf->data[0].second,
                                 buf->size, 0, buf->offset,
                                 sizeof(Point2f));
            }
        }

        ImPlot::EndPlot();
    }
}

}
