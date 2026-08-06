#ifndef CTRLPP_HAS_ARGMIN
#error "the installed backend target reached this translation unit without its compile definition"
#endif

#include <ctrlpp/mpc/argmin_solver.h>

#include <iostream>

int main()
{
    ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp> solver;
    (void)solver;
    std::cout << "ctrlpp integration test PASSED" << std::endl;
    return 0;
}
