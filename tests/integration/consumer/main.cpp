#ifdef TEST_CORE_ONLY
#include <ctrlpp/ctrlpp.h>
#else
#include <ctrlpp/ctrlpp_eigen.h>
#endif

#include <iostream>

int main()
{
#ifdef TEST_CORE_ONLY
    std::cout << "ctrlpp core integration test PASSED\n";
#else
    ctrlpp::DiscreteStateSpace<ctrlpp::EigenLinalgPolicy, double, 2, 1, 1> sys{
        .A = Eigen::Matrix2d::Identity(),
        .B = Eigen::Vector2d::Zero(),
        .C = Eigen::RowVector2d::Zero(),
        .D = Eigen::Matrix<double, 1, 1>::Zero()
    };
    auto x = Eigen::Vector2d::Zero().eval();
    auto u = Eigen::Matrix<double, 1, 1>::Zero().eval();
    auto xn = ctrlpp::propagate(sys, x, u);
    (void)xn;
    std::cout << "ctrlpp eigen integration test PASSED\n";
#endif
    return 0;
}
