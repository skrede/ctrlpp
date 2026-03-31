// so3_exp_log.cpp -- SO(3) exp/log round-trip matching so3_exp_log.m
// Usage: ./so3_exp_log > so3_exp_log_cpp.csv

#include "ctrlpp/lie/so3.h"

#include <Eigen/Geometry>

#include <array>
#include <cstdio>
#include <numbers>

int main()
{
    using Vec3 = Eigen::Vector3d;

    std::array<Vec3, 7> test_vecs = {{
        {0.0, 0.0, 0.0},
        {0.1, 0.0, 0.0},
        {0.0, 0.5, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 1.0, 1.0},
        {std::numbers::pi * 0.99, 0.0, 0.0},
        {0.001, 0.002, 0.003},
    }};

    std::printf("phi_x,phi_y,phi_z,q_w,q_x,q_y,q_z,phi_rt_x,phi_rt_y,phi_rt_z\n");

    for(const auto& phi : test_vecs)
    {
        auto q = ctrlpp::so3::exp(phi);
        auto phi_rt = ctrlpp::so3::log(q);

        std::printf("%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n",
                    phi(0), phi(1), phi(2),
                    q.w(), q.x(), q.y(), q.z(),
                    phi_rt(0), phi_rt(1), phi_rt(2));
    }
}
