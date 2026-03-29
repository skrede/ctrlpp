// cubic_spline_natural.cpp -- Natural cubic spline matching cubic_spline_natural.m
// Usage: ./cubic_spline_natural > cubic_spline_natural_cpp.csv

#include "ctrlpp/trajectory/cubic_spline.h"

#include <cstdio>
#include <vector>

int main()
{
    ctrlpp::cubic_spline<double>::config cfg;
    cfg.times = {0.0, 1.0, 2.5, 4.0, 5.0};
    cfg.positions = {0.0, 1.5, 0.8, 2.0, 1.0};
    cfg.bc = ctrlpp::boundary_condition::natural;

    ctrlpp::cubic_spline<double> spline(cfg);

    constexpr double dt_eval = 0.05;

    std::printf("time,position,velocity,acceleration\n");

    for(double t = cfg.times.front(); t <= cfg.times.back() + dt_eval / 2; t += dt_eval)
    {
        auto pt = spline.evaluate(t);
        std::printf("%.15e,%.15e,%.15e,%.15e\n", t, pt.position(0), pt.velocity(0), pt.acceleration(0));
    }
}
