#include "ctrlpp/implot/app.h"
#include "ctrlpp/implot/sim_harness.h"
#include "ctrlpp/implot/signal_recorder.h"
#include "ctrlpp/implot/scrolling_plot.h"
#include "ctrlpp/implot/static_plot.h"

#include "ctrlpp/ctrlpp.h"

#include "imgui.h"

#include <array>
#include <cmath>
#include <complex>
#include <optional>

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;
    constexpr std::size_t NY = 1;
    constexpr double dt = 0.1;

    using Scalar = double;
    using Mat2 = Eigen::Matrix<Scalar, 2, 2>;
    using Vec2 = Eigen::Matrix<Scalar, 2, 1>;
    using Mat1 = Eigen::Matrix<Scalar, 1, 1>;
    using MatC = Eigen::Matrix<Scalar, 1, 2>;
    using MatB = Eigen::Matrix<Scalar, 2, 1>;
    using MatD = Eigen::Matrix<Scalar, 1, 1>;
    using MatL = Eigen::Matrix<Scalar, 2, 1>;
    using System = ctrlpp::DiscreteStateSpace<Scalar, NX, NU, NY>;

    // Double integrator
    Mat2 A;
    A << 1.0, dt,
         0.0, 1.0;
    MatB B;
    B << 0.5 * dt * dt,
         dt;
    MatC C;
    C << 1.0, 0.0;
    MatD D = MatD::Zero();
    System sys{A, B, C, D};

    // LQR design
    float q00 = 10.0f;
    float q11 = 1.0f;
    float r00 = 1.0f;

    auto make_Q = [](float q0, float q1) {
        Mat2 Q = Mat2::Zero();
        Q(0, 0) = q0;
        Q(1, 1) = q1;
        return Q;
    };
    auto make_R = [](float r) {
        Mat1 R;
        R << static_cast<Scalar>(r);
        return R;
    };

    auto K_opt = ctrlpp::lqr_gain<Scalar, NX, NU>(A, B, make_Q(q00, q11), make_R(r00));
    if (!K_opt) return 1;
    auto lqr = ctrlpp::Lqr<Scalar, NX, NU>(*K_opt);

    // Kalman filter
    float q_proc_diag = 0.01f;
    float r_meas_val = 0.1f;

    auto make_Qproc = [](float q) { return Mat2::Identity() * static_cast<Scalar>(q); };
    auto make_Rmeas = [](float r) {
        Mat1 R;
        R << static_cast<Scalar>(r);
        return R;
    };

    Mat2 P0 = Mat2::Identity();
    Vec2 x0_true;
    x0_true << 1.0, 0.0;
    Vec2 x0_est = Vec2::Zero();

    // Mutable initial conditions for reset
    float x1_true_init = 1.0f;
    float x2_true_init = 0.0f;

    auto kf = ctrlpp::KalmanFilter<Scalar, NX, NU, NY>(
        sys, make_Qproc(q_proc_diag), make_Rmeas(r_meas_val), x0_est, P0);

    // Luenberger observer
    float obs_pole1 = 0.3f;
    float obs_pole2 = 0.4f;

    auto make_poles = [](float p1, float p2) {
        return std::array<std::complex<Scalar>, NX>{
            std::complex<Scalar>{p1, 0.0},
            std::complex<Scalar>{p2, 0.0}
        };
    };

    auto L_opt = ctrlpp::place_observer<Scalar, NX, NY>(A, C, make_poles(obs_pole1, obs_pole2));
    if (!L_opt) return 1;
    auto luenberger = ctrlpp::LuenbergerObserver<Scalar, NX, NU, NY>(sys, *L_opt, x0_est);

    // Simulation state
    Vec2 x_true = x0_true;
    Mat1 u_ctrl = Mat1::Zero();

    ctrlpp::implot::SignalRecorder recorder(4000);

    ctrlpp::implot::SimHarness sim(dt,
        [&](double t) {
            kf.predict(u_ctrl);
            luenberger.predict(u_ctrl);

            x_true = (A * x_true + B * u_ctrl).eval();
            Mat1 z = (C * x_true).eval();

            kf.update(z);
            luenberger.update(z);

            u_ctrl = lqr.compute(kf.state());

            recorder.record("x1_true", t, x_true(0));
            recorder.record("x2_true", t, x_true(1));
            recorder.record("x1_kf", t, kf.state()(0));
            recorder.record("x2_kf", t, kf.state()(1));
            recorder.record("x1_luenberger", t, luenberger.state()(0));
            recorder.record("x2_luenberger", t, luenberger.state()(1));
            recorder.record("control", t, u_ctrl(0, 0));
            recorder.record("kf_err", t, (x_true - kf.state()).norm());
            recorder.record("luenberger_err", t, (x_true - luenberger.state()).norm());
        },
        [&] {
            x_true << static_cast<Scalar>(x1_true_init),
                      static_cast<Scalar>(x2_true_init);
            u_ctrl = Mat1::Zero();
            recorder.clear();

            // Reconstruct observers (no state-reset API)
            kf = ctrlpp::KalmanFilter<Scalar, NX, NU, NY>(
                sys, make_Qproc(q_proc_diag), make_Rmeas(r_meas_val), x0_est, P0);
            luenberger = ctrlpp::LuenbergerObserver<Scalar, NX, NU, NY>(
                sys, L_opt.value_or(MatL::Zero()), x0_est);
        });

    ctrlpp::implot::App app(1400, 800, "LQR + Observer Comparison");

    app.run([&] {
        sim.advance(ImGui::GetIO().DeltaTime);

        // LQR parameters
        ImGui::Begin("LQR Parameters");
        if (ImGui::CollapsingHeader("Simulation", ImGuiTreeNodeFlags_DefaultOpen)) {
            sim.draw_controls();
        }
        ImGui::Separator();
        ImGui::SliderFloat("Q(0,0)", &q00, 0.1f, 100.0f, "%.2f",
            ImGuiSliderFlags_Logarithmic);
        ImGui::SliderFloat("Q(1,1)", &q11, 0.1f, 100.0f, "%.2f",
            ImGuiSliderFlags_Logarithmic);
        ImGui::SliderFloat("R(0,0)", &r00, 0.01f, 10.0f, "%.3f",
            ImGuiSliderFlags_Logarithmic);
        if (ImGui::Button("Recompute LQR Gain")) {
            auto new_K = ctrlpp::lqr_gain<Scalar, NX, NU>(
                A, B, make_Q(q00, q11), make_R(r00));
            if (new_K) {
                K_opt = new_K;
                lqr = ctrlpp::Lqr<Scalar, NX, NU>(*K_opt);
            }
        }
        ImGui::Text("K = [%.4f, %.4f]", K_opt->coeff(0, 0), K_opt->coeff(0, 1));

        ImGui::Separator();
        if (ImGui::Button("Export SVG"))
            ctrlpp::implot::write_svg("lqr_observer.svg", recorder,
                {"x1_true", "x1_kf", "x1_luenberger", "control"});
        ImGui::End();

        // Kalman filter parameters
        ImGui::Begin("Kalman Filter");
        bool kf_changed = false;
        kf_changed |= ImGui::SliderFloat("Q_proc", &q_proc_diag, 0.001f, 1.0f,
            "%.4f", ImGuiSliderFlags_Logarithmic);
        kf_changed |= ImGui::SliderFloat("R_meas", &r_meas_val, 0.01f, 10.0f,
            "%.3f", ImGuiSliderFlags_Logarithmic);
        if (kf_changed)
            kf.set_noise(make_Qproc(q_proc_diag), make_Rmeas(r_meas_val));
        ImGui::Text("Est. error: %.6f", (x_true - kf.state()).norm());
        ImGui::End();

        // Luenberger observer parameters
        ImGui::Begin("Luenberger Observer");
        ImGui::SliderFloat("Pole 1", &obs_pole1, 0.01f, 0.99f);
        ImGui::SliderFloat("Pole 2", &obs_pole2, 0.01f, 0.99f);
        if (ImGui::Button("Recompute Observer Gain")) {
            auto new_L = ctrlpp::place_observer<Scalar, NX, NY>(
                A, C, make_poles(obs_pole1, obs_pole2));
            if (new_L) {
                L_opt = new_L;
                luenberger.set_gain(*L_opt);
            }
        }
        if (L_opt)
            ImGui::Text("L = [%.4f, %.4f]^T", L_opt->coeff(0, 0), L_opt->coeff(1, 0));
        ImGui::Text("Est. error: %.6f", (x_true - luenberger.state()).norm());
        ImGui::End();

        // Initial conditions (take effect on reset)
        ImGui::Begin("Initial Conditions");
        ImGui::SliderFloat("x1_true", &x1_true_init, -5.0f, 5.0f);
        ImGui::SliderFloat("x2_true", &x2_true_init, -5.0f, 5.0f);
        ImGui::TextDisabled("(Applied on Reset)");
        ImGui::End();

        // Plots
        auto t = static_cast<float>(sim.sim_time());
        ctrlpp::implot::scrolling_plot(
            {"Position (x1)", "Time (s)", "x1", 10.0f},
            recorder, {"x1_true", "x1_kf", "x1_luenberger"}, t);
        ctrlpp::implot::scrolling_plot(
            {"Velocity (x2)", "Time (s)", "x2", 10.0f},
            recorder, {"x2_true", "x2_kf", "x2_luenberger"}, t);
        ctrlpp::implot::scrolling_plot(
            {"Control Signal", "Time (s)", "u", 10.0f},
            recorder, {"control"}, t);
        ctrlpp::implot::scrolling_plot(
            {"Estimation Error", "Time (s)", "||e||", 10.0f},
            recorder, {"kf_err", "luenberger_err"}, t);
    });
}
