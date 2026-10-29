#include "ctrlpp/control/care.h"
#include "ctrlpp/control/dare.h"

#include "ctrlpp/detail/riccati_solution.h"
#include "ctrlpp/detail/power_of_two_scale.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <limits>
#include <cstddef>

namespace {

template<int Rows, int Cols>
using matrix = Eigen::Matrix<double, Rows, Cols>;

// fl(s * x) is the exact product iff fma(s, x, -fl(s * x)) is zero: the error-free product transformation
// (Ogita, Rump and Oishi, "Accurate sum and dot product", SIAM J. Sci. Comput. 26(6), 2005, Sec. 3).
template<int Rows, int Cols>
auto is_exact_multiple(const matrix<Rows, Cols> &returned, double scale, const matrix<Rows, Cols> &checked) -> bool
{
    const matrix<Rows, Cols> rounded = scale * checked;
    if(!(returned.array() == rounded.array()).all())
        return false;
    for(int row = 0; row < Rows; ++row)
        for(int col = 0; col < Cols; ++col)
            if(std::fma(scale, checked(row, col), -returned(row, col)) != 0.0)
                return false;
    return true;
}

auto passes_target_definiteness(const matrix<2, 2> &P) -> bool
{
    Eigen::LDLT<matrix<2, 2>> ldlt(P);
    return ldlt.info() == Eigen::Success && ldlt.vectorD().minCoeff() >= ctrlpp::detail::psd_pivot_floor<double, 2>(P);
}

template<std::size_t NX, std::size_t NU>
auto check_care_transport(const matrix<int(NX), int(NX)> &A, const matrix<int(NX), int(NU)> &B, const matrix<int(NX), int(NX)> &Q, const matrix<int(NU), int(NU)> &R) -> void
{
    const double scale                  = ctrlpp::detail::care_weight_scale<double, NX, NU>(Q, R);
    const matrix<int(NX), int(NX)> Q_eq = Q / scale;
    const matrix<int(NU), int(NU)> R_eq = R / scale;
    REQUIRE(ctrlpp::detail::care_weight_scale<double, NX, NU>(Q_eq, R_eq) == 1.0);

    const auto returned = ctrlpp::care<double, NX, NU>(A, B, Q, R);
    const auto checked  = ctrlpp::care<double, NX, NU>(A, B, Q_eq, R_eq);
    REQUIRE(returned.has_value());
    REQUIRE(checked.has_value());
    CHECK(is_exact_multiple(returned->P, scale, checked->P));
}

template<std::size_t NX, std::size_t NU>
auto check_dare_transport(const matrix<int(NX), int(NX)> &A, const matrix<int(NX), int(NU)> &B, const matrix<int(NX), int(NX)> &Q, const matrix<int(NU), int(NU)> &R) -> void
{
    const double scale                  = ctrlpp::detail::dare_weight_scale<double, NX, NU>(Q, R);
    const matrix<int(NX), int(NX)> Q_eq = Q / scale;
    const matrix<int(NU), int(NU)> R_eq = R / scale;
    REQUIRE(ctrlpp::detail::dare_weight_scale<double, NX, NU>(Q_eq, R_eq) == 1.0);

    const auto returned = ctrlpp::dare<double, NX, NU>(A, B, Q, R);
    const auto checked  = ctrlpp::dare<double, NX, NU>(A, B, Q_eq, R_eq);
    REQUIRE(returned.has_value());
    REQUIRE(checked.has_value());
    CHECK(is_exact_multiple(returned->P, scale, checked->P));
    CHECK((returned->K.array() == checked->K.array()).all());
}

constexpr std::array<double, 4> common_scales{3.0, 0.1, 1e-7, 1e150};

// The recorded continuous definiteness inputs of the randomized target, after its clamp and zero-flush:
// A row-major, B, the weight factor row-major, and R.
struct recorded_pose
{
    std::array<double, 4> A;
    std::array<double, 2> B;
    std::array<double, 4> Q_raw;
    double R;
};

constexpr std::array<recorded_pose, 3> recorded_definiteness_inputs{{
        {{-0x1p+1, 0x0p+0, 0x1p+1, 0x0p+0}, {-0x1.fbfbfbf01p-4, -0x1.fbfbfbfbfbfbfp-4}, {-0x1.fbfbfbfbf60ffp-4, -0x1.fbfbfbfbfbfbfp-4, 0x0p+0, 0x0p+0}, 0x1.0624dd2f1a9fcp-8},
        {{0x0p+0, 0x1p+1, 0x1p+1, 0x0p+0}, {0x0p+0, -0x1.fbfbf02p-4}, {-0x1.fbfbfbfbfbfbfp-4, -0x1.fbfbfbfbfbfbfp-4, -0x1.fbfbfbf002fbfp-4, -0x1.fbfbfbfp-4}, 0x1.0624dd2f1a9fcp-8},
        {{0x0p+0, 0x1p+1, 0x0p+0, -0x1p+1}, {0x0p+0, -0x1p+1}, {-0x1.fbfbfbfe3e3e3p-4, -0x1.fbfbfbfbfbfbfp-4, -0x1p+1, -0x1p+1}, 0x1.0624dd2f1a9fcp-8},
}};

}

TEST_CASE("The power-of-two ceiling leaves a power of two unchanged", "[riccati][scale]")
{
    using limits = std::numeric_limits<double>;
    CHECK(ctrlpp::detail::power_of_two_ceiling(1.0) == 1.0);
    CHECK(ctrlpp::detail::power_of_two_ceiling(0.5) == 0.5);
    CHECK(ctrlpp::detail::power_of_two_ceiling(limits::denorm_min()) == limits::denorm_min());
    CHECK(ctrlpp::detail::power_of_two_ceiling(std::ldexp(1.0, limits::max_exponent - 1)) == std::ldexp(1.0, limits::max_exponent - 1));
}

TEST_CASE("The power-of-two ceiling rounds any other magnitude up to the next power of two", "[riccati][scale]")
{
    CHECK(ctrlpp::detail::power_of_two_ceiling(3.0) == 4.0);
    CHECK(ctrlpp::detail::power_of_two_ceiling(0.75) == 1.0);
    CHECK(ctrlpp::detail::power_of_two_ceiling(0.1) == std::ldexp(1.0, -3));
    CHECK(ctrlpp::detail::power_of_two_ceiling(1.0 + std::numeric_limits<double>::epsilon()) == 2.0);
}

TEST_CASE("The power-of-two ceiling is clamped to the largest finite power of two", "[riccati][scale]")
{
    using limits             = std::numeric_limits<double>;
    const double largest_two = std::ldexp(1.0, limits::max_exponent - 1);
    CHECK(ctrlpp::detail::power_of_two_ceiling(1e308) == largest_two);
    CHECK(ctrlpp::detail::power_of_two_ceiling(limits::max()) == largest_two);
}

TEST_CASE("The power-of-two ceiling keeps the zero, NaN and infinite results the solvers decline on", "[riccati][scale]")
{
    using limits = std::numeric_limits<double>;
    CHECK(ctrlpp::detail::power_of_two_ceiling(0.0) == 1.0);
    CHECK(ctrlpp::detail::power_of_two_ceiling(limits::quiet_NaN()) == 1.0);
    CHECK(ctrlpp::detail::power_of_two_ceiling(limits::infinity()) == limits::infinity());
}

TEST_CASE("The power-of-two ceiling covers the single-precision range", "[riccati][scale]")
{
    using limits            = std::numeric_limits<float>;
    const float largest_two = std::ldexp(1.0f, limits::max_exponent - 1);
    CHECK(ctrlpp::detail::power_of_two_ceiling(limits::denorm_min()) == limits::denorm_min());
    CHECK(ctrlpp::detail::power_of_two_ceiling(largest_two) == largest_two);
    CHECK(ctrlpp::detail::power_of_two_ceiling(limits::max()) == largest_two);
}

TEST_CASE("The continuous solution is returned as an exact multiple of the equilibrated one", "[riccati][scale]")
{
    matrix<2, 2> A;
    A << 0.0, 1.0, 0.0, 0.0;
    const matrix<2, 1> B{0.0, 1.0};
    for(const double c : common_scales)
    {
        CAPTURE(c);
        check_care_transport<2, 1>(A, B, c * matrix<2, 2>::Identity(), matrix<1, 1>::Constant(c));
    }

    matrix<2, 2> Q;
    Q << 3.0, 1.0, 1.0, 2.0;
    check_care_transport<2, 1>(A, B, Q, matrix<1, 1>::Constant(5.0));
}

TEST_CASE("The discrete solution is returned as an exact multiple of the equilibrated one with the same gain", "[riccati][scale]")
{
    const matrix<1, 1> A = matrix<1, 1>::Constant(0.5);
    const matrix<1, 1> B = matrix<1, 1>::Constant(1.0);
    for(const double c : common_scales)
    {
        CAPTURE(c);
        check_dare_transport<1, 1>(A, B, matrix<1, 1>::Constant(c), matrix<1, 1>::Constant(c));
    }

    matrix<2, 2> A2;
    A2 << 1.0, 1.0, 0.0, 1.0;
    const matrix<2, 1> B2{0.5, 1.0};
    matrix<2, 2> Q;
    Q << 3.0, 1.0, 1.0, 2.0;
    check_dare_transport<2, 1>(A2, B2, Q, matrix<1, 1>::Constant(5.0));
}

TEST_CASE("The recorded continuous definiteness inputs reach the same verdict at both scales", "[riccati][scale]")
{
    for(const recorded_pose &pose : recorded_definiteness_inputs)
    {
        const matrix<2, 2> A     = Eigen::Map<const Eigen::Matrix<double, 2, 2, Eigen::RowMajor>>(pose.A.data());
        const matrix<2, 1> B     = Eigen::Map<const matrix<2, 1>>(pose.B.data());
        const matrix<2, 2> Q_raw = Eigen::Map<const Eigen::Matrix<double, 2, 2, Eigen::RowMajor>>(pose.Q_raw.data());
        const matrix<2, 2> Q     = Q_raw.transpose() * Q_raw;
        const matrix<1, 1> R     = matrix<1, 1>::Constant(pose.R);

        const double scale  = ctrlpp::detail::care_weight_scale<double, 2, 1>(Q, R);
        const auto returned = ctrlpp::care<double, 2, 1>(A, B, Q, R);
        const auto checked  = ctrlpp::care<double, 2, 1>(A, B, matrix<2, 2>(Q / scale), matrix<1, 1>(R / scale));
        REQUIRE(returned.has_value() == checked.has_value());
        if(returned.has_value())
        {
            CHECK(is_exact_multiple(returned->P, scale, checked->P));
            CHECK(passes_target_definiteness(returned->P));
        }
        else
        {
            CHECK(returned.error() == checked.error());
        }
    }
}
