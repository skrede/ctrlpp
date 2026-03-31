#ifndef HPP_GUARD_TESTS_UNIT_HARDENING_HELPERS_H
#define HPP_GUARD_TESTS_UNIT_HARDENING_HELPERS_H

#include "ctrlpp/types.h"
#include "ctrlpp/model/state_space.h"

#include <cmath>
#include <limits>

namespace ctrlpp::test
{

template <typename Scalar, std::size_t N>
auto nan_vector() -> Vector<Scalar, N>
{
    return Vector<Scalar, N>::Constant(std::numeric_limits<Scalar>::quiet_NaN());
}

template <typename Scalar, std::size_t N>
auto inf_vector() -> Vector<Scalar, N>
{
    return Vector<Scalar, N>::Constant(std::numeric_limits<Scalar>::infinity());
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
auto nan_matrix() -> Matrix<Scalar, Rows, Cols>
{
    return Matrix<Scalar, Rows, Cols>::Constant(std::numeric_limits<Scalar>::quiet_NaN());
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
auto inf_matrix() -> Matrix<Scalar, Rows, Cols>
{
    return Matrix<Scalar, Rows, Cols>::Constant(std::numeric_limits<Scalar>::infinity());
}

template <typename Scalar>
auto ill_conditioned_2x2(Scalar condition_number) -> Matrix<Scalar, 2, 2>
{
    auto m = Matrix<Scalar, 2, 2>::Zero().eval();
    m(0, 0) = Scalar{1};
    m(1, 1) = Scalar{1} / condition_number;
    return m;
}

template <typename Scalar, std::size_t N>
auto ill_conditioned_nxn(Scalar condition_number) -> Matrix<Scalar, N, N>
{
    auto m = Matrix<Scalar, N, N>::Zero().eval();
    for(std::size_t i = 0; i < N; ++i)
        m(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(i)) =
            std::pow(condition_number, -static_cast<Scalar>(i) / static_cast<Scalar>(N - 1));
    return m;
}

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
auto near_singular_system(Scalar epsilon) -> discrete_state_space<Scalar, NX, NU, NY>
{
    auto A = Matrix<Scalar, NX, NX>::Identity().eval();
    A(0, 0) = epsilon;

    auto B = Matrix<Scalar, NX, NU>::Zero().eval();
    for(std::size_t i = 0; i < std::min(NX, NU); ++i)
        B(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(i)) = Scalar{1};

    auto C = Matrix<Scalar, NY, NX>::Zero().eval();
    for(std::size_t i = 0; i < std::min(NY, NX); ++i)
        C(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(i)) = Scalar{1};

    auto D = Matrix<Scalar, NY, NU>::Zero();

    return {A, B, C, D};
}

template <typename Scalar, std::size_t N>
auto zero_vector() -> Vector<Scalar, N>
{
    return Vector<Scalar, N>::Zero();
}

template <typename Scalar, std::size_t Rows, std::size_t Cols>
auto zero_matrix() -> Matrix<Scalar, Rows, Cols>
{
    return Matrix<Scalar, Rows, Cols>::Zero();
}

}

#endif
