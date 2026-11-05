#ifndef HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_PF_STORAGE_H
#define HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_PF_STORAGE_H

#include "ctrlpp/estimation/detail/filter_configuration.h"

#include "ctrlpp/estimation/particle_filter_types.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/detail/covariance_ops.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <random>
#include <cstddef>
#include <numbers>
#include <utility>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t NP, typename Dynamics, typename Measurement, typename Resampler, typename Rng,
        typename Likelihood>
class pf_storage
{
public:
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int ny = static_cast<int>(NY);

    using state_vector_t = Vector<Scalar, NX>;
    using input_vector_t = Vector<Scalar, NU>;
    using output_vector_t = Vector<Scalar, NY>;

protected:
    // Factors S with S S^T equal to Q (process) and to P0 (prior).
    struct sampling_factors
    {
        Matrix<Scalar, NX, NX> process;
        Matrix<Scalar, NX, NX> prior;
    };

    Dynamics m_dynamics;
    Measurement m_measurement;
    Resampler m_resampler;
    Rng m_rng;
    Likelihood m_likelihood;

    Matrix<Scalar, NX, NX> m_Q;
    Matrix<Scalar, NX, NX> m_Q_sqrt;
    Matrix<Scalar, NY, NY> m_R;
    Matrix<Scalar, NY, NY> m_R_inv;
    Scalar m_log_det_2piR;
    Scalar m_ess_threshold;
    Scalar m_roughening_scale;
    extraction_method m_extraction;
    weight_representation m_weight_mode;
    pf_health m_health;

    std::array<state_vector_t, NP> m_particles;
    // The set before the current resampling: the resampling gathers from it,
    // and a roughening that leaves a coordinate non-finite is undone from it.
    std::array<state_vector_t, NP> m_previous_particles;
    std::array<Scalar, NP> m_log_weights;
    std::array<Scalar, NP> m_linear_weights;
    // Staged so a refused update leaves the weights untouched; a member rather
    // than a local keeps the update's stack frame independent of NP.
    std::array<Scalar, NP> m_staged_weights;

    pf_storage(Dynamics&& dynamics, Measurement&& measurement, pf_config<Scalar, NX, NU, NY>&& config, sampling_factors&& factors, Resampler&& resampler, Rng&& rng,
            Likelihood&& likelihood)
        : m_dynamics{std::move(dynamics)}
        , m_measurement{std::move(measurement)}
        , m_resampler{std::move(resampler)}
        , m_rng{std::move(rng)}
        , m_likelihood{std::move(likelihood)}
        , m_Q{std::move(config.Q)}
        , m_Q_sqrt{std::move(factors.process)}
        , m_R{std::move(config.R)}
        , m_R_inv{m_R.colPivHouseholderQr().inverse()}
        , m_log_det_2piR{compute_log_det_2piR()}
        , m_ess_threshold{config.ess_threshold < Scalar{0} ? static_cast<Scalar>(NP) / Scalar{2} : config.ess_threshold}
        , m_roughening_scale{config.roughening_scale}
        , m_extraction{config.extraction}
        , m_weight_mode{config.weights}
        , m_health{pf_health::ok}
        , m_previous_particles{}
        , m_log_weights{}
        , m_linear_weights{}
        , m_staged_weights{}
    {
        initialize_particles(config.x0, factors.prior);
    }

    Scalar compute_log_det_2piR() const
    {
        Scalar log_det = m_R.colPivHouseholderQr().logAbsDeterminant();
        return static_cast<Scalar>(NY) * std::log(Scalar{2} * std::numbers::pi_v<Scalar>) + log_det;
    }

    void initialize_particles(const state_vector_t& x0, const Matrix<Scalar, NX, NX>& P0_sqrt)
    {
        std::normal_distribution<Scalar> normal(Scalar{0}, Scalar{1});

        for(std::size_t i = 0; i < NP; ++i)
        {
            state_vector_t noise;
            for(std::size_t d = 0; d < NX; ++d)
                noise(static_cast<Eigen::Index>(d)) = normal(m_rng);
            m_particles[i] = x0 + P0_sqrt * noise;
        }

        Scalar log_uniform = -std::log(static_cast<Scalar>(NP));
        Scalar lin_uniform = Scalar{1} / static_cast<Scalar>(NP);
        m_log_weights.fill(log_uniform);
        m_linear_weights.fill(lin_uniform);
    }

    static ctrlpp::expected<sampling_factors, filter_error> validated_sampling_factors(pf_config<Scalar, NX, NU, NY>& config)
    {
        const auto Q = detail::checked_covariance(config.Q, detail::process_noise_faults);
        if(!Q)
            return ctrlpp::unexpected(Q.error());
        const auto R = detail::checked_covariance(config.R, detail::measurement_noise_faults);
        if(!R)
            return ctrlpp::unexpected(R.error());
        if(!detail::nonsingular_spectrum(*R))
            return ctrlpp::unexpected(filter_error::singular_measurement_noise);
        if(!config.x0.allFinite())
            return ctrlpp::unexpected(filter_error::non_finite_initial_state);
        const auto P0 = detail::checked_covariance(config.P0, detail::initial_covariance_faults);
        if(!P0)
            return ctrlpp::unexpected(P0.error());
        return sampling_factors{sampling_factor(config.Q, *Q), sampling_factor(config.P0, *P0)};
    }

    // Eigen's LLT leaves the input's entries in the columns after a non-positive pivot, so its
    // factor is read only after a successful factorization.
    static Matrix<Scalar, NX, NX> sampling_factor(const Matrix<Scalar, NX, NX>& M, const detail::psd_decomposition_result<Scalar, nx>& decomposition)
    {
        Eigen::LLT<Matrix<Scalar, NX, NX>> llt(M);
        if(llt.info() == Eigen::Success)
            return llt.matrixL();
        return detail::psd_factor(decomposition);
    }

    state_vector_t sample_process_noise()
    {
        std::normal_distribution<Scalar> normal(Scalar{0}, Scalar{1});
        state_vector_t noise;
        for(std::size_t d = 0; d < NX; ++d)
            noise(static_cast<Eigen::Index>(d)) = normal(m_rng);
        return m_Q_sqrt * noise;
    }

    void propagate_particles(const input_vector_t& u)
    {
        for(std::size_t i = 0; i < NP; ++i)
        {
            m_particles[i] = m_dynamics(m_particles[i], u);
            m_particles[i] += sample_process_noise();
        }
    }
};

}

#endif
