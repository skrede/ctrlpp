#ifndef HPP_GUARD_CTRLPP_EIGEN_LINALG_H
#define HPP_GUARD_CTRLPP_EIGEN_LINALG_H

#include <ctrlpp/linalg_policy.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <array>
#include <limits>
#include <complex>
#include <cstddef>
#include <utility>

namespace ctrlpp {

struct EigenLinalgPolicy {
    template<typename Scalar, std::size_t Rows, std::size_t Cols>
    using matrix_type = Eigen::Matrix<Scalar, static_cast<int>(Rows), static_cast<int>(Cols)>;

    template<typename Scalar, std::size_t Rows>
    using vector_type = Eigen::Matrix<Scalar, static_cast<int>(Rows), 1>;

    template<typename Scalar, int RA, int CA, int CB>
    static auto multiply(const Eigen::Matrix<Scalar, RA, CA>& A,
                          const Eigen::Matrix<Scalar, CA, CB>& B)
        -> Eigen::Matrix<Scalar, RA, CB>
    {
        return (A * B).eval();
    }

    template<typename Scalar, int Rows, int Cols>
    static auto multiply(const Eigen::Matrix<Scalar, Rows, Cols>& A,
                          const Eigen::Matrix<Scalar, Cols, 1>& v)
        -> Eigen::Matrix<Scalar, Rows, 1>
    {
        return (A * v).eval();
    }

    template<typename Scalar, int Rows, int Cols>
    static auto add(const Eigen::Matrix<Scalar, Rows, Cols>& A,
                     const Eigen::Matrix<Scalar, Rows, Cols>& B)
        -> Eigen::Matrix<Scalar, Rows, Cols>
    {
        return (A + B).eval();
    }

    template<typename Scalar, int Rows, int Cols>
    static auto subtract(const Eigen::Matrix<Scalar, Rows, Cols>& A,
                          const Eigen::Matrix<Scalar, Rows, Cols>& B)
        -> Eigen::Matrix<Scalar, Rows, Cols>
    {
        return (A - B).eval();
    }

    template<typename Scalar, int N>
    static auto solve(const Eigen::Matrix<Scalar, N, N>& A,
                       const Eigen::Matrix<Scalar, N, 1>& b)
        -> Eigen::Matrix<Scalar, N, 1>
    {
        return A.colPivHouseholderQr().solve(b).eval();
    }

    template<typename Scalar, int N, int Cols>
    static auto solve(const Eigen::Matrix<Scalar, N, N>& A,
                       const Eigen::Matrix<Scalar, N, Cols>& B)
        -> Eigen::Matrix<Scalar, N, Cols>
    {
        return A.colPivHouseholderQr().solve(B).eval();
    }

    template<typename Scalar, int Rows, int Cols>
    static auto transpose(const Eigen::Matrix<Scalar, Rows, Cols>& A)
        -> Eigen::Matrix<Scalar, Cols, Rows>
    {
        return A.transpose().eval();
    }

    template<typename Scalar, std::size_t N>
    static auto identity() -> matrix_type<Scalar, N, N>
    {
        return matrix_type<Scalar, N, N>::Identity();
    }

    template<typename Scalar, int Rows, int Cols>
    static auto get_element(const Eigen::Matrix<Scalar, Rows, Cols>& A,
                            std::size_t row, std::size_t col) -> Scalar
    {
        return A(static_cast<int>(row), static_cast<int>(col));
    }

    template<typename Scalar, int Rows, int Cols>
    static void set_element(Eigen::Matrix<Scalar, Rows, Cols>& A,
                            std::size_t row, std::size_t col, Scalar val)
    {
        A(static_cast<int>(row), static_cast<int>(col)) = val;
    }

    template<typename Scalar, std::size_t N>
    static auto zeros() -> matrix_type<Scalar, N, N>
    {
        return matrix_type<Scalar, N, N>::Zero();
    }

    template<typename Scalar, int N>
    static auto trace(const Eigen::Matrix<Scalar, N, N>& A) -> Scalar
    {
        return A.trace();
    }

    template<typename Scalar, int N>
    static auto eigenvalues(const Eigen::Matrix<Scalar, N, N>& A)
        -> std::array<std::complex<Scalar>, static_cast<std::size_t>(N)>
    {
        Eigen::EigenSolver<Eigen::Matrix<Scalar, N, N>> solver(A, false);
        auto evals = solver.eigenvalues();
        std::array<std::complex<Scalar>, static_cast<std::size_t>(N)> result;
        for (int i = 0; i < N; ++i)
            result[static_cast<std::size_t>(i)] = evals(i);
        return result;
    }

    template<typename Scalar, int N>
    static auto schur(const Eigen::Matrix<Scalar, N, N>& A)
        -> std::pair<Eigen::Matrix<Scalar, N, N>, Eigen::Matrix<Scalar, N, N>>
    {
        Eigen::RealSchur<Eigen::Matrix<Scalar, N, N>> decomposition(A);
        return {decomposition.matrixT().eval(), decomposition.matrixU().eval()};
    }

    template<typename Scalar, int Rows, int Cols>
    static auto rank(const Eigen::Matrix<Scalar, Rows, Cols>& A) -> std::size_t
    {
        return static_cast<std::size_t>(A.colPivHouseholderQr().rank());
    }
};

static_assert(
    std::numeric_limits<int>::max() > 0,
    "Eigen dimension overflow guard: ensure int can represent dimensions. "
    "EigenLinalgPolicy uses static_cast<int> for Eigen dimension parameters.");

namespace detail {

template<typename Scalar, std::size_t N, typename Policy>
concept SchurPolicy = LinalgPolicy<Policy> && requires(
    typename Policy::template matrix_type<Scalar, N, N> A
) {
    { Policy::schur(A) } -> std::same_as<
        std::pair<typename Policy::template matrix_type<Scalar, N, N>,
                  typename Policy::template matrix_type<Scalar, N, N>>>;
};

template<typename Scalar, std::size_t R, std::size_t C, typename Policy>
concept RankPolicy = LinalgPolicy<Policy> && requires(
    typename Policy::template matrix_type<Scalar, R, C> A
) {
    { Policy::rank(A) } -> std::convertible_to<std::size_t>;
};

}

}

#endif
