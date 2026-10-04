#include "driver.h"

int main(int argc, char** argv)
{
    return ctrlpp::bench::qp_solvers::run_driver(ctrlpp::bench::qp_solvers::cohort::sparse, argc, argv);
}
