#include "ctrlpp/control/pid.h"

int main()
{
#ifdef CTRLPP_CANARY_VELOCITY_ANTI_WINDUP
    using controller_t = ctrlpp::pid<double, 1, ctrlpp::velocity_form, ctrlpp::anti_windup<CTRLPP_CANARY_VELOCITY_ANTI_WINDUP>>;
#else
    using controller_t = ctrlpp::pid<double, 1, ctrlpp::velocity_form>;
#endif
    const controller_t::config_type cfg{};
    auto controller = controller_t::create(cfg);
#ifdef CTRLPP_CANARY_VELOCITY_SET_INTEGRAL
    return controller && controller->set_integral(controller_t::vector_t::Zero()) ? 0 : 1;
#else
    return controller.has_value() ? 0 : 1;
#endif
}
