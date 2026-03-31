// butterworth_filter.cpp -- 4th-order Butterworth matching butterworth_filter.m
// Usage: ./butterworth_filter > butterworth_filter_cpp.csv

#include "ctrlpp/dsp/biquad.h"

#include <cstdio>

int main()
{
    auto filt = ctrlpp::make_butterworth<4>(50.0, 1000.0);

    constexpr int n_samples = 200;

    std::printf("sample,response\n");

    for(int i = 0; i < n_samples; ++i)
    {
        double output = filt.process(1.0);
        std::printf("%.15e,%.15e\n", static_cast<double>(i), output);
    }
}
