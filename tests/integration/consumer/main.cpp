#include <ctrlpp/ctrlpp.h>

#include <iostream>

int main()
{
    ctrlpp::DiscreteStateSpace<double, 2, 1, 1> sys{
        .A = Eigen::Matrix2d::Identity(),
        .B = Eigen::Vector2d::Zero(),
        .C = Eigen::RowVector2d::Zero(),
        .D = Eigen::Matrix<double, 1, 1>::Zero()
    };
    auto x = Eigen::Vector2d::Zero().eval();
    auto u = Eigen::Matrix<double, 1, 1>::Zero().eval();
    auto xn = ctrlpp::propagate(sys, x, u);
    (void)xn;
    std::cout << "ctrlpp integration test PASSED" << std::endl;
    return 0;
}
