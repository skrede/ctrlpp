#include <ctrlpp/model/propagate.h>
#include <ctrlpp/model/state_space.h>

#include <iostream>

int main()
{
    ctrlpp::discrete_state_space<double, 2, 1, 1> sys{Eigen::Matrix2d::Identity(), Eigen::Vector2d::Zero(), Eigen::RowVector2d::Zero(), Eigen::Matrix<double, 1, 1>::Zero()};
    auto x = Eigen::Vector2d::Zero().eval();
    auto u = Eigen::Matrix<double, 1, 1>::Zero().eval();
    auto xn = ctrlpp::propagate(sys, x, u);
    (void)xn;
    std::cout << "ctrlpp integration test PASSED" << std::endl;
    return 0;
}
