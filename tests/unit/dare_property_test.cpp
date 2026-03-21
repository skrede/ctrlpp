#include "ctrlpp/dare.h"
#include "ctrlpp/lqr.h"
#include "ctrlpp/analysis.h"

#include <catch2/catch_test_macros.hpp>
#include <rapidcheck.h>
#include <rapidcheck/catch.h>

#include <cmath>
#include <complex>
#include <cstddef>

namespace {

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr int nx = static_cast<int>(NX);
constexpr int nu = static_cast<int>(NU);

using Mat2 = Eigen::Matrix<double, nx, nx>;
using Mat21 = Eigen::Matrix<double, nx, nu>;
using Mat1 = Eigen::Matrix<double, nu, nu>;

auto bounded_double(double lo, double hi) -> rc::Gen<double>
{
    return rc::gen::map(rc::gen::inRange(0, 1000000), [lo, hi](int x) {
        return lo + (hi - lo) * (static_cast<double>(x) / 1000000.0);
    });
}

auto gen_random_matrix(int rows, int cols, double lo, double hi) -> rc::Gen<Eigen::MatrixXd>
{
    return rc::gen::map(
        rc::gen::container<std::vector<double>>(
            static_cast<std::size_t>(rows * cols),
            bounded_double(lo, hi)),
        [rows, cols](const std::vector<double>& vals) {
            Eigen::MatrixXd m(rows, cols);
            for (int i = 0; i < rows * cols; ++i)
                m(i / cols, i % cols) = vals[static_cast<std::size_t>(i)];
            return m;
        });
}

auto gen_psd_matrix() -> rc::Gen<Mat2>
{
    return rc::gen::map(gen_random_matrix(nx, nx, -5.0, 5.0), [](const Eigen::MatrixXd& Ld) {
        Mat2 L = Ld.cast<double>();
        return Mat2{(L.transpose() * L + 1e-6 * Mat2::Identity()).eval()};
    });
}

auto gen_pd_matrix_1x1() -> rc::Gen<Mat1>
{
    return rc::gen::map(bounded_double(0.1, 10.0), [](double v) {
        Mat1 m;
        m << v;
        return m;
    });
}

struct ControllableSystem {
    Mat2 A;
    Mat21 B;
};

auto gen_controllable_system() -> rc::Gen<ControllableSystem>
{
    return rc::gen::suchThat(
        rc::gen::apply(
            [](const Eigen::MatrixXd& Ad, const Eigen::MatrixXd& Bd) {
                ControllableSystem sys;
                sys.A = Ad.cast<double>();
                sys.B = Bd.cast<double>();
                return sys;
            },
            gen_random_matrix(nx, nx, -2.0, 2.0),
            gen_random_matrix(nx, nu, -2.0, 2.0)),
        [](const ControllableSystem& sys) {
            return ctrlpp::is_controllable<double, NX, NU>(sys.A, sys.B);
        });
}

}

TEST_CASE("dare riccati identity holds for random controllable systems", "[dare][property]")
{
    rc::prop("riccati residual is near zero for valid solutions", [] {
        auto sys = *gen_controllable_system();
        auto Q = *gen_psd_matrix();
        auto R = *gen_pd_matrix_1x1();

        auto P_opt = ctrlpp::dare<double, NX, NU>(sys.A, sys.B, Q, R);
        if (!P_opt) {
            RC_SUCCEED("dare returned nullopt for this system");
        }

        auto& P = *P_opt;
        auto BtP = (sys.B.transpose() * P).eval();
        auto S = (R + BtP * sys.B).eval();
        auto S_inv_BtPA = S.colPivHouseholderQr().solve(BtP * sys.A).eval();

        // Riccati residual: A^T P A - P - A^T P B (R + B^T P B)^{-1} B^T P A + Q
        Mat2 residual = sys.A.transpose() * P * sys.A - P
                      - sys.A.transpose() * P * sys.B * S_inv_BtPA + Q;

        double rel_tol = 1e-8 * std::max(P.norm(), 1.0);
        RC_ASSERT(residual.norm() < rel_tol);
    });
}

TEST_CASE("dare solution is symmetric positive semi-definite", "[dare][property]")
{
    rc::prop("P is symmetric and PSD when dare succeeds", [] {
        auto sys = *gen_controllable_system();
        auto Q = *gen_psd_matrix();
        auto R = *gen_pd_matrix_1x1();

        auto P_opt = ctrlpp::dare<double, NX, NU>(sys.A, sys.B, Q, R);
        if (!P_opt) {
            RC_SUCCEED("dare returned nullopt");
        }

        auto& P = *P_opt;

        // Symmetry
        RC_ASSERT((P - P.transpose()).norm() < 1e-12);

        // PSD: all eigenvalues >= -epsilon
        Eigen::SelfAdjointEigenSolver<Mat2> eig(P);
        for (int i = 0; i < nx; ++i) {
            RC_ASSERT(eig.eigenvalues()(i) >= -1e-10);
        }
    });
}

TEST_CASE("lqr closed-loop eigenvalues inside unit circle", "[lqr][property]")
{
    rc::prop("closed-loop system is stable", [] {
        auto sys = *gen_controllable_system();
        auto Q = *gen_psd_matrix();
        auto R = *gen_pd_matrix_1x1();

        auto K_opt = ctrlpp::lqr_gain<double, NX, NU>(sys.A, sys.B, Q, R);
        if (!K_opt) {
            RC_SUCCEED("lqr_gain returned nullopt");
        }

        auto& K = *K_opt;
        Mat2 Acl = sys.A - sys.B * K;

        Eigen::EigenSolver<Mat2> eig(Acl, false);
        for (int i = 0; i < nx; ++i) {
            RC_ASSERT(std::abs(eig.eigenvalues()(i)) < 1.0);
        }
    });
}

TEST_CASE("dare robustness - degenerate inputs", "[dare][property]")
{
    rc::prop("dare never crashes on degenerate inputs", [] {
        // Generate matrices with potentially extreme values
        auto A_gen = *gen_random_matrix(nx, nx, -100.0, 100.0);
        auto B_gen = *gen_random_matrix(nx, nu, -100.0, 100.0);
        auto Q_gen = *gen_psd_matrix();
        auto R_gen = *gen_pd_matrix_1x1();

        Mat2 A = A_gen.cast<double>();
        Mat21 B = B_gen.cast<double>();

        auto P_opt = ctrlpp::dare<double, NX, NU>(A, B, Q_gen, R_gen);

        if (P_opt) {
            auto& P = *P_opt;
            // If returned, all entries must be finite
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < nx; ++j)
                    RC_ASSERT(std::isfinite(P(i, j)));
        }
        // nullopt is acceptable -- the point is no crash
    });
}
