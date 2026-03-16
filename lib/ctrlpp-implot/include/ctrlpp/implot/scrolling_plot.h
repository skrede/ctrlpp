#ifndef HPP_GUARD_CTRLPP_IMPLOT_SCROLLING_PLOT_H
#define HPP_GUARD_CTRLPP_IMPLOT_SCROLLING_PLOT_H

#include "ctrlpp/implot/signal_recorder.h"

#include <initializer_list>
#include <string_view>

namespace ctrlpp::implot {

struct PlotConfig {
    const char* title{"Signal"};
    const char* x_label{"Time (s)"};
    const char* y_label{"Value"};
    float history{10.0f};
    float height{0.0f};
};

void scrolling_plot(const PlotConfig& config, const SignalRecorder& recorder,
                    std::initializer_list<const char*> signal_names,
                    float current_time);

}

#endif
