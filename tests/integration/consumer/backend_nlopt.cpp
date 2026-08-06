#ifndef CTRLPP_HAS_NLOPT
#error "the installed backend target reached this translation unit without its compile definition"
#endif

#include <ctrlpp/mpc/nlopt_solver.h>

#include <iostream>

int main()
{
    ctrlpp::nlopt_solver<double> solver;
    (void)solver;
    std::cout << "ctrlpp integration test PASSED" << std::endl;
    return 0;
}
