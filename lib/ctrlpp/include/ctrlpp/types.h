#ifndef HPP_GUARD_CTRLPP_TYPES_H
#define HPP_GUARD_CTRLPP_TYPES_H

#include <Eigen/Dense>
#include <cstddef>

namespace ctrlpp {

template<typename Scalar, std::size_t Rows, std::size_t Cols>
using Matrix = Eigen::Matrix<Scalar, static_cast<int>(Rows), static_cast<int>(Cols)>;

template<typename Scalar, std::size_t Rows>
using Vector = Eigen::Matrix<Scalar, static_cast<int>(Rows), 1>;

}

#endif
