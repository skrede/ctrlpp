#ifndef HPP_GUARD_TESTS_UNIT_HARDENING_HELPERS_H
#define HPP_GUARD_TESTS_UNIT_HARDENING_HELPERS_H

#include "ctrlpp/types.h"

#include "ctrlpp/model/state_space.h"

#include "ctrlpp/detail/riccati_solution.h"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <utility>
#include <algorithm>
#include <type_traits>

namespace ctrlpp::test
{

/// @brief Assert that a fallible controller cycle produced a command, and hand
/// the command back.
///
/// This is the presence assertion, not an unchecked unwrap: a rejected cycle
/// fails the enclosing test case at the point of the call instead of reading
/// through an empty result. Tests that mean to OBSERVE a rejection assert on
/// the result directly and never come through here, so the helper cannot hide
/// one.
template <typename Result>
auto commanded(const Result& result) -> std::remove_cvref_t<decltype(*result)>
{
    REQUIRE(result.has_value());
    return *result;
}

/// @brief Assert that a fallible construction produced an object, and hand the
/// object back.
///
/// The same presence assertion as `commanded`, applied to the construction
/// channel: a rejected configuration fails the enclosing test case where it was
/// configured, rather than being unwrapped through an empty result or papered
/// over with a substituted object. Cases that mean to OBSERVE a rejection assert
/// on the result directly and never come through here, so the helper cannot hide
/// one.
template <typename Result>
auto constructed(Result&& result) -> std::remove_cvref_t<decltype(*result)>
{
    REQUIRE(result.has_value());
    return *std::forward<Result>(result);
}

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

/// @brief The gain the discrete Riccati solution implies:
/// K = (R + B'PB)^{-1} B'PA.
template <typename Scalar, std::size_t NX, std::size_t NU>
auto riccati_gain(const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NX, NU>& B,
                  const Matrix<Scalar, NU, NU>& R, const Matrix<Scalar, NX, NX>& P)
    -> Matrix<Scalar, NU, NX>
{
    auto const BtP = (B.transpose() * P).eval();
    auto const S = (R + BtP * B).eval();
    return S.colPivHouseholderQr().solve(BtP * A).eval();
}

/// @brief Residual of the discrete algebraic Riccati equation at P, paired with
/// the largest of the four terms whose cancellation produces it.
///
/// The scale is NOT the solution norm alone. The residual is a sum of four
/// matrices that cancel almost exactly, and its rounding is set by the largest
/// of them: for a state matrix of magnitude much above or below one, A'PA
/// dwarfs P and a budget written against P alone would be measuring the wrong
/// quantity.
template <typename Scalar, std::size_t NX, std::size_t NU>
struct riccati_residual_result
{
    Scalar norm{};
    Scalar scale{};
};

template <typename Scalar, std::size_t NX, std::size_t NU>
auto riccati_residual(const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NX, NU>& B,
                      const Matrix<Scalar, NX, NX>& Q, const Matrix<Scalar, NU, NU>& R,
                      const Matrix<Scalar, NX, NX>& P) -> riccati_residual_result<Scalar, NX, NU>
{
    auto const K = riccati_gain<Scalar, NX, NU>(A, B, R, P);
    auto const AtPA = (A.transpose() * P * A).eval();
    auto const AtPBK = (A.transpose() * P * B * K).eval();
    auto const residual = (AtPA - P - AtPBK + Q).eval();

    Scalar const scale = std::max({AtPA.norm(), P.norm(), AtPBK.norm(), Q.norm()});
    return {residual.norm(), scale};
}

/// @brief The rounded-operation count along the longest chain producing one
/// entry of the Riccati residual.
///
/// A QUOTE, not a copy. The enumeration lives with the library's own bound in
/// `ctrlpp/detail/riccati_solution.h`; this name exists only so the test call
/// sites read in test vocabulary. The hand-copy that used to sit here was a
/// second definition of a bound the library also defines, which is exactly the
/// drift path the shared positive-semi-definiteness floor was consolidated to
/// close.
template <std::size_t NX, std::size_t NU>
constexpr int riccati_residual_ops = ctrlpp::detail::dare_residual_ops<NX, NU>;

/// @brief Whether a solved Riccati pose retains more than half the scalar
/// type's significand, decided by the library's own rule.
///
/// Also a quote. The estimator, the margin and the comparison are all defined
/// once, at `ctrlpp::detail::riccati_forward_error_verdict`; an anchor that
/// re-spelled any of the three would be asserting against its own copy of the
/// contract rather than against the contract.
template <typename Scalar, std::size_t NX, std::size_t NU>
auto riccati_accuracy_of(const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NX, NU>& B,
                         const Matrix<Scalar, NU, NU>& R, const Matrix<Scalar, NX, NX>& Q,
                         const Matrix<Scalar, NX, NX>& P) -> ctrlpp::detail::riccati_accuracy
{
    auto const K = riccati_gain<Scalar, NX, NU>(A, B, R, P);
    auto const closed_loop = (A - B * K).eval();
    auto const AtPA = (A.transpose() * P * A).eval();
    auto const AtPBK = (A.transpose() * P * B * K).eval();
    auto const residual = (AtPA - P - AtPBK + Q).eval();
    return ctrlpp::detail::riccati_forward_error_verdict<Scalar, int(NX)>(closed_loop, residual, P);
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
