#ifndef CTRLPP_HAS_OSQP
#error "the installed backend target reached this translation unit without its compile definition"
#endif

#include <ctrlpp/mpc/osqp_solver.h>

#include <iostream>

int main()
{
    ctrlpp::osqp_solver solver{ctrlpp::qp_preset::speed};
    (void)solver;
    std::cout << "ctrlpp integration test PASSED" << std::endl;
    return 0;
}
