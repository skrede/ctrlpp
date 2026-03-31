// fir_filter.cpp -- FIR filter response matching fir_filter.m
// Usage: ./fir_filter > fir_filter_cpp.csv

#include "ctrlpp/dsp/fir.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <numbers>

int main()
{
    ctrlpp::fir<double, 5> filt(std::array<double, 5>{0.2, 0.2, 0.2, 0.2, 0.2});

    constexpr int n_samples = 100;

    std::printf("sample,input,output\n");

    for(int i = 0; i < n_samples; ++i)
    {
        double x = 0.0;
        if(i >= 9)
            x = 1.0;
        x += 0.3 * std::sin(2.0 * std::numbers::pi * 0.05 * i);

        double y = filt.process(x);
        std::printf("%.15e,%.15e,%.15e\n", static_cast<double>(i), x, y);
    }
}
