#include "perf_common.h"

#include "ctrlpp/mpc/argmin_solver.h"

#include <vector>
#include <cstdlib>
#include <iostream>
#include <exception>

int main(int argc, char** argv)
{
    using namespace ctrlpp::argmin_perf;

    run_options opts;
    try
    {
        opts = parse_options(argc, argv);
    }
    catch(const std::exception& e)
    {
        std::cerr << "perf_argmin_slsqp: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    using Solver = ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp>;

    std::vector<trace_row> trace;
    auto* trace_sink = opts.trace_path.empty() ? nullptr : &trace;

    auto [x_final, total_cost] = run_pendulum_closed_loop<Solver>(
        opts, trace_sink, "argmin_slsqp");

    if(trace_sink)
        dump_trace(opts.trace_path, trace);

    std::cout << "x_final=[" << x_final(0) << ',' << x_final(1) << ']'
              << " total_cost=" << total_cost << std::endl;

    return EXIT_SUCCESS;
}
