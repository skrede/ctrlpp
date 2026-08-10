#include "ctrlpp/sysid/recursive_arx.h"

int main()
{
#ifdef CTRLPP_CANARY_FIVE_PARAMETER_ARITY
    const auto estimator = ctrlpp::recursive_arx<double, 2, 1, 1, 1>::create();
#else
    const auto estimator = ctrlpp::recursive_arx<double, 2, 1>::create();
#endif
    return estimator ? 0 : 1;
}
