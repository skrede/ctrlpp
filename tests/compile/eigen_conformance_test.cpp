#include "ctrlpp/eigen_linalg.h"

// EigenLinalgPolicy satisfies LinalgPolicy (concept checks at double)
static_assert(ctrlpp::LinalgPolicy<ctrlpp::EigenLinalgPolicy>);

// Verify type aliases resolve correctly for float
static_assert(std::is_same_v<
    ctrlpp::EigenLinalgPolicy::matrix_type<float, 3, 3>,
    Eigen::Matrix<float, 3, 3>>);

static_assert(std::is_same_v<
    ctrlpp::EigenLinalgPolicy::vector_type<float, 3>,
    Eigen::Matrix<float, 3, 1>>);

// Verify type aliases resolve correctly for double
static_assert(std::is_same_v<
    ctrlpp::EigenLinalgPolicy::matrix_type<double, 3, 3>,
    Eigen::Matrix<double, 3, 3>>);

static_assert(std::is_same_v<
    ctrlpp::EigenLinalgPolicy::vector_type<double, 3>,
    Eigen::Matrix<double, 3, 1>>);

int main() { return 0; }
