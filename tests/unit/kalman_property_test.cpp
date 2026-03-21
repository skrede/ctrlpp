#include "ctrlpp/kalman.h"
#include "ctrlpp/analysis.h"
#include "ctrlpp/state_space.h"

#include <catch2/catch_test_macros.hpp>
#include <rapidcheck.h>
#include <rapidcheck/catch.h>

#include <cmath>
#include <cstddef>

namespace {

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr std::size_t NY = 1;
constexpr int nx = static_cast<int>(NX);
constexpr int nu = static_cast<int>(NU);
constexpr int ny = static_cast<int>(NY);

using Mat2 = Eigen::Matrix<double, nx, nx>;
using Mat21 = Eigen::Matrix<double, nx, nu>;
using Mat12 = Eigen::Matrix<double, ny, nx>;
using Mat11 = Eigen::Matrix<double, ny, ny>;
using Mat11u = Eigen::Matrix<double, nu, nu>;
using Vec2 = Eigen::Matrix<double, nx, 1>;
using Vec1 = Eigen::Matrix<double, ny, 1>;
using Vec1u = Eigen::Matrix<double, nu, 1>;

auto bounded_double(double lo, double hi) -> rc::Gen<double>
{
    return rc::gen::map(rc::gen::inRange(0, 1000000), [lo, hi](int x) {
        return lo + (hi - lo) * (static_cast<double>(x) / 1000000.0);
    });
}

struct MassSpringDamper {
    Mat2 A;
    Mat21 B;
    Mat12 C;
    Mat11 D;
};

auto gen_mass_spring_damper() -> rc::Gen<MassSpringDamper>
{
    return rc::gen::suchThat(
        rc::gen::apply(
            [](double m, double k, double c, double dt) {
                // Continuous: A = [[0,1],[-k/m,-c/m]], B=[[0],[1/m]]
                // Forward Euler discretization: Ad = I + A*dt, Bd = B*dt
                Mat2 Ac;
                Ac << 0.0, 1.0,
                      -k / m, -c / m;
                Mat21 Bc;
                Bc << 0.0, 1.0 / m;

                MassSpringDamper sys;
                sys.A = Mat2::Identity() + Ac * dt;
                sys.B = Bc * dt;
                sys.C << 1.0, 0.0;
                sys.D << 0.0;
                return sys;
            },
            bounded_double(0.5, 10.0),   // mass
            bounded_double(0.1, 50.0),   // spring
            bounded_double(0.1, 10.0),   // damper
            bounded_double(0.001, 0.05)  // dt (small for Euler stability)
        ),
        [](const MassSpringDamper& sys) {
            return ctrlpp::is_observable<double, NX, NY>(sys.A, sys.C);
        });
}

auto gen_psd_2x2() -> rc::Gen<Mat2>
{
    return rc::gen::apply(
        [](double a, double b, double c, double d) {
            Mat2 L;
            L << a, b, c, d;
            return Mat2{(L.transpose() * L + 1e-6 * Mat2::Identity()).eval()};
        },
        bounded_double(-3.0, 3.0),
        bounded_double(-3.0, 3.0),
        bounded_double(-3.0, 3.0),
        bounded_double(-3.0, 3.0));
}

auto gen_pd_1x1() -> rc::Gen<Mat11>
{
    return rc::gen::map(bounded_double(0.01, 10.0), [](double v) {
        Mat11 m;
        m << v;
        return m;
    });
}

}

TEST_CASE("kalman covariance remains psd throughout filtering", "[kalman][property]")
{
    rc::prop("covariance eigenvalues non-negative at every step", [] {
        auto sys = *gen_mass_spring_damper();
        auto Q = *gen_psd_2x2();
        auto R = *gen_pd_1x1();

        ctrlpp::discrete_state_space<double, NX, NU, NY> dss{sys.A, sys.B, sys.C, sys.D};
        Vec2 x0 = Vec2::Zero();
        Mat2 P0 = Mat2::Identity();

        ctrlpp::kalman_filter<double, NX, NU, NY> kf(dss, {.Q = Q, .R = R, .x0 = x0, .P0 = P0});

        Vec2 x_true = Vec2::Zero();
        Vec1u u = Vec1u::Zero();

        for (int step = 0; step < 50; ++step) {
            kf.predict(u);
            x_true = (sys.A * x_true + sys.B * u).eval();

            Vec1 z = sys.C * x_true;
            kf.update(z);

            auto& P = kf.covariance();
            Eigen::SelfAdjointEigenSolver<Mat2> eig(P);
            for (int i = 0; i < nx; ++i) {
                RC_ASSERT(eig.eigenvalues()(i) >= -1e-10);
            }
        }
    });
}

TEST_CASE("kalman covariance is symmetric", "[kalman][property]")
{
    rc::prop("P - P^T near zero at every step", [] {
        auto sys = *gen_mass_spring_damper();
        auto Q = *gen_psd_2x2();
        auto R = *gen_pd_1x1();

        ctrlpp::discrete_state_space<double, NX, NU, NY> dss{sys.A, sys.B, sys.C, sys.D};
        Vec2 x0 = Vec2::Zero();
        Mat2 P0 = Mat2::Identity();

        ctrlpp::kalman_filter<double, NX, NU, NY> kf(dss, {.Q = Q, .R = R, .x0 = x0, .P0 = P0});

        Vec2 x_true = Vec2::Zero();
        Vec1u u = Vec1u::Zero();

        for (int step = 0; step < 50; ++step) {
            kf.predict(u);
            x_true = (sys.A * x_true + sys.B * u).eval();

            Vec1 z = sys.C * x_true;
            kf.update(z);

            auto& P = kf.covariance();
            RC_ASSERT((P - P.transpose()).norm() < 1e-12);
        }
    });
}

TEST_CASE("kalman innovation bounded for stable system", "[kalman][property]")
{
    rc::prop("innovation does not diverge", [] {
        auto sys = *gen_mass_spring_damper();
        auto Q = *gen_psd_2x2();
        auto R = *gen_pd_1x1();

        ctrlpp::discrete_state_space<double, NX, NU, NY> dss{sys.A, sys.B, sys.C, sys.D};
        Vec2 x0 = Vec2::Zero();
        Mat2 P0 = 10.0 * Mat2::Identity();

        ctrlpp::kalman_filter<double, NX, NU, NY> kf(dss, {.Q = Q, .R = R, .x0 = x0, .P0 = P0});

        // Offset initial true state to create observable innovation
        Vec2 x_true;
        x_true << 1.0, 0.0;
        Vec1u u = Vec1u::Zero();

        double max_innov_first_half = 0.0;
        double max_innov_second_half = 0.0;
        constexpr int total_steps = 100;

        for (int step = 0; step < total_steps; ++step) {
            kf.predict(u);
            x_true = (sys.A * x_true + sys.B * u).eval();

            Vec1 z = sys.C * x_true;
            kf.update(z);

            double innov_norm = kf.innovation().norm();
            if (step < total_steps / 2)
                max_innov_first_half = std::max(max_innov_first_half, innov_norm);
            else
                max_innov_second_half = std::max(max_innov_second_half, innov_norm);
        }

        // Second half innovation should not exceed first half (convergence)
        RC_ASSERT(max_innov_second_half <= max_innov_first_half + 1e-6);
    });
}

TEST_CASE("kalman robustness - extreme noise parameters", "[kalman][property]")
{
    rc::prop("no NaN/Inf with extreme Q/R ratios", [] {
        auto sys = *gen_mass_spring_damper();

        auto q_scale = *bounded_double(1e-10, 1e6);
        auto r_scale = *bounded_double(1e-10, 1e6);

        Mat2 Q = q_scale * Mat2::Identity();
        Mat11 R;
        R << r_scale;

        ctrlpp::discrete_state_space<double, NX, NU, NY> dss{sys.A, sys.B, sys.C, sys.D};
        Vec2 x0 = Vec2::Zero();
        Mat2 P0 = Mat2::Identity();

        ctrlpp::kalman_filter<double, NX, NU, NY> kf(dss, {.Q = Q, .R = R, .x0 = x0, .P0 = P0});

        Vec2 x_true = Vec2::Zero();
        Vec1u u = Vec1u::Zero();

        for (int step = 0; step < 20; ++step) {
            kf.predict(u);
            x_true = (sys.A * x_true + sys.B * u).eval();

            Vec1 z = sys.C * x_true;
            kf.update(z);

            auto& state = kf.state();
            auto& P = kf.covariance();

            for (int i = 0; i < nx; ++i) {
                RC_ASSERT(std::isfinite(state(i)));
                for (int j = 0; j < nx; ++j)
                    RC_ASSERT(std::isfinite(P(i, j)));
            }
        }
    });
}
