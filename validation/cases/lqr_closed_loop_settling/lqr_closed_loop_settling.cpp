#include "control_loop_demo.h"

#include <cstdio>
#include <cstddef>

int main()
{
    auto demo = ctrlpp::control_loop_demo<double>::make();
    if(!demo.has_value())
    {
        std::fprintf(stderr, "the gain design declined the embedded plant: %s\n", ctrlpp::describe(demo.error()));
        return 1;
    }

    for(std::size_t k = 0; k < ctrlpp::kGoldenSteps; ++k)
        demo->step();

    std::printf("K_00,K_01,x0_final,x1_final,final_norm\n");
    std::printf("%.15e,%.15e,%.15e,%.15e,%.15e\n", demo->K(0, 0), demo->K(0, 1), demo->x(0), demo->x(1), demo->x.norm());
    return 0;
}
