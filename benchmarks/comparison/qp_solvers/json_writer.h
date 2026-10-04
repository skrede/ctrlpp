#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_JSON_WRITER_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_JSON_WRITER_H

#include <cmath>
#include <format>
#include <string>
#include <vector>
#include <cstdint>
#include <ostream>
#include <string_view>

namespace ctrlpp::bench::qp_solvers
{

/// A streaming writer for the handful of JSON objects the drivers hand to the run script. Floats are
/// written in shortest round-trip form and always as floats, so 1.0 stays a float for a reader that
/// compares JSON types.
class json_writer
{
public:
    explicit json_writer(std::ostream& out)
        : out_(out)
        , first_{}
    {
    }

    void open()
    {
        separate();
        out_ << '{';
        first_.push_back(true);
    }

    void open(std::string_view key)
    {
        name(key);
        out_ << '{';
        first_.push_back(true);
    }

    void close()
    {
        out_ << '}';
        first_.pop_back();
    }

    void text(std::string_view key, std::string_view value)
    {
        name(key);
        out_ << '"' << value << '"';
    }

    void integer(std::string_view key, std::int64_t value)
    {
        name(key);
        out_ << value;
    }

    void flag(std::string_view key, bool value)
    {
        name(key);
        out_ << (value ? "true" : "false");
    }

    void number(std::string_view key, double value)
    {
        name(key);
        out_ << float_text(value);
    }

    static auto float_text(double value) -> std::string
    {
        if(!std::isfinite(value))
            return "null";
        std::string text = std::format("{}", value);
        if(text.find_first_of(".e") == std::string::npos)
            text += ".0";
        return text;
    }

private:
    std::ostream& out_;
    std::vector<bool> first_;

    void separate()
    {
        if(!first_.empty() && !first_.back())
            out_ << ',';
        if(!first_.empty())
            first_.back() = false;
    }

    void name(std::string_view key)
    {
        separate();
        out_ << '"' << key << "\":";
    }
};

}

#endif
