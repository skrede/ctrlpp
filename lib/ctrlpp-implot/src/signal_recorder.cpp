#include "ctrlpp/implot/signal_recorder.h"

#include <algorithm>

namespace ctrlpp::implot {

RingBuffer::RingBuffer(int capacity)
    : max_size{capacity}
{
    data.resize(static_cast<std::size_t>(capacity));
}

void RingBuffer::push(float t, float v)
{
    if (size < max_size) {
        data[static_cast<std::size_t>(size)] = ImVec2(t, v);
        ++size;
        offset = size % max_size;
    } else {
        data[static_cast<std::size_t>(offset)] = ImVec2(t, v);
        offset = (offset + 1) % max_size;
    }
}

void RingBuffer::clear()
{
    size = 0;
    offset = 0;
}

SignalRecorder::SignalRecorder(int capacity)
    : capacity_{capacity}
{
}

void SignalRecorder::record(std::string_view name, double time, double value)
{
    auto [it, _] = channels_.try_emplace(std::string(name), capacity_);
    it->second.push(static_cast<float>(time), static_cast<float>(value));
}

const RingBuffer* SignalRecorder::channel(std::string_view name) const
{
    auto it = channels_.find(std::string(name));
    if (it != channels_.end()) {
        return &it->second;
    }
    return nullptr;
}

void SignalRecorder::clear()
{
    for (auto& [_, buf] : channels_) {
        buf.clear();
    }
}

void SignalRecorder::clear(std::string_view name)
{
    auto it = channels_.find(std::string(name));
    if (it != channels_.end()) {
        it->second.clear();
    }
}

std::vector<std::string> SignalRecorder::channel_names() const
{
    std::vector<std::string> names;
    names.reserve(channels_.size());
    for (const auto& [name, _] : channels_) {
        names.push_back(name);
    }
    return names;
}

}
