#ifndef HPP_GUARD_CTRLPP_IMPLOT_SIGNAL_RECORDER_H
#define HPP_GUARD_CTRLPP_IMPLOT_SIGNAL_RECORDER_H

#include <cstddef>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ctrlpp::implot {

using Point2f = std::pair<float, float>;

struct RingBuffer {
    std::vector<Point2f> data;
    int size{0};
    int offset{0};
    int max_size;

    explicit RingBuffer(int capacity);

    void push(float t, float v);
    void clear();
};

class SignalRecorder {
public:
    explicit SignalRecorder(int capacity = 2000);

    void record(std::string_view name, double time, double value);
    [[nodiscard]] const RingBuffer* channel(std::string_view name) const;
    void clear();
    void clear(std::string_view name);
    [[nodiscard]] std::vector<std::string> channel_names() const;

private:
    int capacity_;
    std::unordered_map<std::string, RingBuffer> channels_;
};

}

#endif
