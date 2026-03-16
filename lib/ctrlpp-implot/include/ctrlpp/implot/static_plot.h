#ifndef HPP_GUARD_CTRLPP_IMPLOT_STATIC_PLOT_H
#define HPP_GUARD_CTRLPP_IMPLOT_STATIC_PLOT_H

#include "ctrlpp/implot/signal_recorder.h"

#include <initializer_list>

namespace ctrlpp::implot {

struct SvgConfig {
    const char* title{"Signal"};
    const char* x_label{"Time (s)"};
    const char* y_label{"Value"};
    int width{800};
    int height{400};
};

void write_svg(const char* path,
               const SignalRecorder& recorder,
               std::initializer_list<const char*> signals,
               const SvgConfig& config = {});

}

#endif
