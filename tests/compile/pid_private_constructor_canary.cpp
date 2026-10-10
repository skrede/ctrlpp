#include "ctrlpp/control/pid.h"

int main()
{
    const ctrlpp::pid<double, 1>::config_type cfg{};
#ifdef CTRLPP_CANARY_PID_DIRECT_CONSTRUCTION
    ctrlpp::pid<double, 1> controller(cfg);
    return controller.integral()(0) == 0.0 ? 0 : 1;
#else
    const auto controller = ctrlpp::pid<double, 1>::create(cfg);
    return controller.has_value() ? 0 : 1;
#endif
}
