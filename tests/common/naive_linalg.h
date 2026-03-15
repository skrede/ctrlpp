#ifndef HPP_GUARD_CPPCTRL_TEST_NAIVE_LINALG_H
#define HPP_GUARD_CPPCTRL_TEST_NAIVE_LINALG_H

#include <array>
#include <cstddef>

struct NaiveLinalg {
    template<typename Scalar, std::size_t Rows, std::size_t Cols>
    using matrix_type = std::array<std::array<Scalar, Cols>, Rows>;

    template<typename Scalar, std::size_t Rows>
    using vector_type = std::array<Scalar, Rows>;

    template<typename Scalar, std::size_t RA, std::size_t CA, std::size_t CB>
    static constexpr auto multiply(const matrix_type<Scalar, RA, CA>& A,
                                   const matrix_type<Scalar, CA, CB>& B)
        -> matrix_type<Scalar, RA, CB>
    {
        matrix_type<Scalar, RA, CB> result{};
        for (std::size_t i = 0; i < RA; ++i)
            for (std::size_t j = 0; j < CB; ++j)
                for (std::size_t k = 0; k < CA; ++k)
                    result[i][j] += A[i][k] * B[k][j];
        return result;
    }

    template<typename Scalar, std::size_t Rows, std::size_t Cols>
    static constexpr auto multiply(const matrix_type<Scalar, Rows, Cols>& A,
                                   const vector_type<Scalar, Cols>& v)
        -> vector_type<Scalar, Rows>
    {
        vector_type<Scalar, Rows> result{};
        for (std::size_t i = 0; i < Rows; ++i)
            for (std::size_t k = 0; k < Cols; ++k)
                result[i] += A[i][k] * v[k];
        return result;
    }

    template<typename Scalar, std::size_t N>
    static constexpr auto solve(const matrix_type<Scalar, N, N>& A,
                                const vector_type<Scalar, N>& b)
        -> vector_type<Scalar, N>
    {
        if constexpr (N == 2) {
            auto det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
            return vector_type<Scalar, 2>{
                (b[0] * A[1][1] - b[1] * A[0][1]) / det,
                (A[0][0] * b[1] - A[1][0] * b[0]) / det
            };
        } else {
            return vector_type<Scalar, N>{};
        }
    }

    template<typename Scalar, std::size_t Rows, std::size_t Cols>
    static constexpr auto transpose(const matrix_type<Scalar, Rows, Cols>& A)
        -> matrix_type<Scalar, Cols, Rows>
    {
        matrix_type<Scalar, Cols, Rows> result{};
        for (std::size_t i = 0; i < Rows; ++i)
            for (std::size_t j = 0; j < Cols; ++j)
                result[j][i] = A[i][j];
        return result;
    }

    template<typename Scalar, std::size_t N>
    static constexpr auto identity() -> matrix_type<Scalar, N, N>
    {
        matrix_type<Scalar, N, N> result{};
        for (std::size_t i = 0; i < N; ++i)
            result[i][i] = Scalar{1};
        return result;
    }

    template<typename Scalar, std::size_t Rows, std::size_t Cols>
    static constexpr auto add(const matrix_type<Scalar, Rows, Cols>& A,
                              const matrix_type<Scalar, Rows, Cols>& B)
        -> matrix_type<Scalar, Rows, Cols>
    {
        matrix_type<Scalar, Rows, Cols> result{};
        for (std::size_t i = 0; i < Rows; ++i)
            for (std::size_t j = 0; j < Cols; ++j)
                result[i][j] = A[i][j] + B[i][j];
        return result;
    }

    template<typename Scalar, std::size_t Rows>
    static constexpr auto add(const vector_type<Scalar, Rows>& a,
                              const vector_type<Scalar, Rows>& b)
        -> vector_type<Scalar, Rows>
    {
        vector_type<Scalar, Rows> result{};
        for (std::size_t i = 0; i < Rows; ++i)
            result[i] = a[i] + b[i];
        return result;
    }

    template<typename Scalar, std::size_t Rows, std::size_t Cols>
    static constexpr auto subtract(const matrix_type<Scalar, Rows, Cols>& A,
                                   const matrix_type<Scalar, Rows, Cols>& B)
        -> matrix_type<Scalar, Rows, Cols>
    {
        matrix_type<Scalar, Rows, Cols> result{};
        for (std::size_t i = 0; i < Rows; ++i)
            for (std::size_t j = 0; j < Cols; ++j)
                result[i][j] = A[i][j] - B[i][j];
        return result;
    }

    template<typename Scalar, std::size_t Rows>
    static constexpr auto subtract(const vector_type<Scalar, Rows>& a,
                                   const vector_type<Scalar, Rows>& b)
        -> vector_type<Scalar, Rows>
    {
        vector_type<Scalar, Rows> result{};
        for (std::size_t i = 0; i < Rows; ++i)
            result[i] = a[i] - b[i];
        return result;
    }

    template<typename Scalar, std::size_t N, std::size_t Cols>
    static constexpr auto solve(const matrix_type<Scalar, N, N>& A,
                                const matrix_type<Scalar, N, Cols>& B)
        -> matrix_type<Scalar, N, Cols>
    {
        if constexpr (N == 2) {
            auto det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
            matrix_type<Scalar, N, Cols> result{};
            for (std::size_t j = 0; j < Cols; ++j) {
                result[0][j] = (B[0][j] * A[1][1] - B[1][j] * A[0][1]) / det;
                result[1][j] = (A[0][0] * B[1][j] - A[1][0] * B[0][j]) / det;
            }
            return result;
        } else {
            return matrix_type<Scalar, N, Cols>{};
        }
    }
};

#endif
