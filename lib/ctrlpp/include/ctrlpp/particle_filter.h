#ifndef HPP_GUARD_CTRLPP_PARTICLE_FILTER_H
#define HPP_GUARD_CTRLPP_PARTICLE_FILTER_H

#include "ctrlpp/types.h"
#include "ctrlpp/observer_policy.h"
#include "ctrlpp/mpc/dynamics_model.h"
#include "ctrlpp/mpc/measurement_model.h"
#include "ctrlpp/resampling/resampling_strategy.h"
#include "ctrlpp/resampling/systematic_resampling.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <random>
#include <type_traits>
#include <utility>

namespace ctrlpp {

enum class extraction_method {
    weighted_mean,
    map
};

enum class weight_representation {
    log,
    linear
};

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct pf_config {
    Matrix<Scalar, NX, NX> Q{Matrix<Scalar, NX, NX>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Vector<Scalar, NX> x0{Vector<Scalar, NX>::Zero()};
    Matrix<Scalar, NX, NX> P0{Matrix<Scalar, NX, NX>::Identity()};
    Scalar ess_threshold{Scalar{-1}};
    Scalar roughening_scale{Scalar{0.2}};
    extraction_method extraction{extraction_method::weighted_mean};
    weight_representation weights{weight_representation::log};
};

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY,
         std::size_t NP, typename Dynamics, typename Measurement,
         typename Resampler = systematic_resampling,
         typename Rng = std::mt19937_64>
requires dynamics_model<Dynamics, Scalar, NX, NU> &&
         measurement_model<Measurement, Scalar, NX, NY> &&
         std::uniform_random_bit_generator<Rng> &&
         resampling_strategy<Resampler, Rng, NP>
class particle_filter {
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int ny = static_cast<int>(NY);

public:
    using observer_tag = struct pf_tag;
    using state_vector_t = Vector<Scalar, NX>;
    using input_vector_t = Vector<Scalar, NU>;
    using output_vector_t = Vector<Scalar, NY>;

    particle_filter(Dynamics dynamics, Measurement measurement,
                    pf_config<Scalar, NX, NU, NY> config, Rng rng = Rng{})
        : dynamics_{std::move(dynamics)}
        , measurement_{std::move(measurement)}
        , resampler_{}
        , rng_{std::move(rng)}
        , Q_{std::move(config.Q)}
        , R_{std::move(config.R)}
        , R_inv_{R_.colPivHouseholderQr().inverse()}
        , log_det_2piR_{compute_log_det_2piR()}
        , ess_threshold_{config.ess_threshold < Scalar{0}
                         ? static_cast<Scalar>(NP) / Scalar{2}
                         : config.ess_threshold}
        , roughening_scale_{config.roughening_scale}
        , extraction_{config.extraction}
        , weight_mode_{config.weights}
    {
        initialize_particles(config.x0, config.P0);
    }

    particle_filter(Dynamics dynamics, Measurement measurement,
                    pf_config<Scalar, NX, NU, NY> config,
                    Resampler resampler, Rng rng = Rng{})
        : dynamics_{std::move(dynamics)}
        , measurement_{std::move(measurement)}
        , resampler_{std::move(resampler)}
        , rng_{std::move(rng)}
        , Q_{std::move(config.Q)}
        , R_{std::move(config.R)}
        , R_inv_{R_.colPivHouseholderQr().inverse()}
        , log_det_2piR_{compute_log_det_2piR()}
        , ess_threshold_{config.ess_threshold < Scalar{0}
                         ? static_cast<Scalar>(NP) / Scalar{2}
                         : config.ess_threshold}
        , roughening_scale_{config.roughening_scale}
        , extraction_{config.extraction}
        , weight_mode_{config.weights}
    {
        initialize_particles(config.x0, config.P0);
    }

    void predict(const input_vector_t& u)
    {
        for (std::size_t i = 0; i < NP; ++i) {
            particles_[i] = dynamics_(particles_[i], u);
            particles_[i] += sample_process_noise();
        }
    }

    void update(const output_vector_t& z)
    {
        if (weight_mode_ == weight_representation::log) {
            update_log(z);
        } else {
            update_linear(z);
        }
    }

    [[nodiscard]] auto state() const -> state_vector_t
    {
        if (extraction_ == extraction_method::map) {
            return map_estimate();
        }
        return weighted_mean();
    }

    [[nodiscard]] auto weighted_mean() const -> state_vector_t
    {
        state_vector_t mean = state_vector_t::Zero();
        if (weight_mode_ == weight_representation::log) {
            auto lin = log_to_linear();
            for (std::size_t i = 0; i < NP; ++i) {
                mean += lin[i] * particles_[i];
            }
        } else {
            for (std::size_t i = 0; i < NP; ++i) {
                mean += linear_weights_[i] * particles_[i];
            }
        }
        return mean;
    }

    [[nodiscard]] auto map_estimate() const -> state_vector_t
    {
        std::size_t best = 0;
        if (weight_mode_ == weight_representation::log) {
            for (std::size_t i = 1; i < NP; ++i) {
                if (log_weights_[i] > log_weights_[best]) {
                    best = i;
                }
            }
        } else {
            for (std::size_t i = 1; i < NP; ++i) {
                if (linear_weights_[i] > linear_weights_[best]) {
                    best = i;
                }
            }
        }
        return particles_[best];
    }

    [[nodiscard]] auto particles() const -> const std::array<state_vector_t, NP>& { return particles_; }

private:
    Dynamics dynamics_;
    Measurement measurement_;
    Resampler resampler_;
    Rng rng_;

    Matrix<Scalar, NX, NX> Q_;
    Matrix<Scalar, NY, NY> R_;
    Matrix<Scalar, NY, NY> R_inv_;
    Scalar log_det_2piR_;
    Scalar ess_threshold_;
    Scalar roughening_scale_;
    extraction_method extraction_;
    weight_representation weight_mode_;

    std::array<state_vector_t, NP> particles_;
    std::array<Scalar, NP> log_weights_{};
    std::array<Scalar, NP> linear_weights_{};

    // Cholesky factor of Q for process noise sampling
    Eigen::Matrix<Scalar, nx, nx> Q_L_{};

    [[nodiscard]] auto compute_log_det_2piR() const -> Scalar
    {
        Scalar log_det = R_.colPivHouseholderQr().logAbsDeterminant();
        return static_cast<Scalar>(NY) * std::log(Scalar{2} * std::numbers::pi_v<Scalar>) + log_det;
    }

    void initialize_particles(const state_vector_t& x0, const Matrix<Scalar, NX, NX>& P0)
    {
        // LLT decomposition of P0 for sampling
        Eigen::LLT<Eigen::Matrix<Scalar, nx, nx>> llt_P0(P0);
        Eigen::Matrix<Scalar, nx, nx> sqrt_P0 = llt_P0.matrixL();

        // LLT decomposition of Q for process noise
        Eigen::LLT<Eigen::Matrix<Scalar, nx, nx>> llt_Q(Q_);
        Q_L_ = llt_Q.matrixL();

        std::normal_distribution<Scalar> normal(Scalar{0}, Scalar{1});

        for (std::size_t i = 0; i < NP; ++i) {
            state_vector_t noise;
            for (std::size_t d = 0; d < NX; ++d) {
                noise(static_cast<int>(d)) = normal(rng_);
            }
            particles_[i] = x0 + sqrt_P0 * noise;
        }

        // Uniform weights
        Scalar log_uniform = -std::log(static_cast<Scalar>(NP));
        Scalar lin_uniform = Scalar{1} / static_cast<Scalar>(NP);
        log_weights_.fill(log_uniform);
        linear_weights_.fill(lin_uniform);
    }

    [[nodiscard]] auto sample_process_noise() -> state_vector_t
    {
        std::normal_distribution<Scalar> normal(Scalar{0}, Scalar{1});
        state_vector_t noise;
        for (std::size_t d = 0; d < NX; ++d) {
            noise(static_cast<int>(d)) = normal(rng_);
        }
        return Q_L_ * noise;
    }

    [[nodiscard]] auto log_likelihood(const output_vector_t& z,
                                      const output_vector_t& z_pred) const -> Scalar
    {
        output_vector_t innov = z - z_pred;
        Scalar mahal = (innov.transpose() * R_inv_ * innov)(0, 0);
        return Scalar{-0.5} * mahal - Scalar{0.5} * log_det_2piR_;
    }

    void update_log(const output_vector_t& z)
    {
        // Update log weights
        for (std::size_t i = 0; i < NP; ++i) {
            output_vector_t z_pred = measurement_(particles_[i]);
            log_weights_[i] += log_likelihood(z, z_pred);
        }

        // Log-sum-exp normalization
        Scalar max_log_w = *std::max_element(log_weights_.begin(), log_weights_.end());
        Scalar sum_exp = Scalar{0};
        for (std::size_t i = 0; i < NP; ++i) {
            sum_exp += std::exp(log_weights_[i] - max_log_w);
        }
        Scalar log_sum = max_log_w + std::log(sum_exp);
        for (auto& lw : log_weights_) {
            lw -= log_sum;
        }

        // Compute ESS
        Scalar sum_w2 = Scalar{0};
        for (std::size_t i = 0; i < NP; ++i) {
            Scalar w = std::exp(log_weights_[i]);
            sum_w2 += w * w;
        }
        Scalar ess = Scalar{1} / sum_w2;

        if (ess < ess_threshold_) {
            // Convert to linear for resampling
            auto lin = log_to_linear();

            std::array<std::size_t, NP> indices{};
            resampler_.resample(lin, indices, rng_);

            reindex_particles(indices);
            apply_roughening();

            // Reset to uniform
            Scalar log_uniform = -std::log(static_cast<Scalar>(NP));
            log_weights_.fill(log_uniform);
        }
    }

    void update_linear(const output_vector_t& z)
    {
        for (std::size_t i = 0; i < NP; ++i) {
            output_vector_t z_pred = measurement_(particles_[i]);
            Scalar ll = log_likelihood(z, z_pred);
            linear_weights_[i] *= std::exp(ll);
        }

        // Normalize
        Scalar sum = Scalar{0};
        for (auto w : linear_weights_) sum += w;
        if (sum > Scalar{0}) {
            for (auto& w : linear_weights_) w /= sum;
        }

        // Compute ESS
        Scalar sum_w2 = Scalar{0};
        for (auto w : linear_weights_) sum_w2 += w * w;
        Scalar ess = Scalar{1} / sum_w2;

        if (ess < ess_threshold_) {
            std::array<std::size_t, NP> indices{};
            resampler_.resample(linear_weights_, indices, rng_);

            reindex_particles(indices);
            apply_roughening();

            Scalar uniform = Scalar{1} / static_cast<Scalar>(NP);
            linear_weights_.fill(uniform);
        }
    }

    [[nodiscard]] auto log_to_linear() const -> std::array<Scalar, NP>
    {
        std::array<Scalar, NP> lin;
        for (std::size_t i = 0; i < NP; ++i) {
            lin[i] = std::exp(log_weights_[i]);
        }
        return lin;
    }

    void reindex_particles(const std::array<std::size_t, NP>& indices)
    {
        std::array<state_vector_t, NP> resampled;
        for (std::size_t i = 0; i < NP; ++i) {
            resampled[i] = particles_[indices[i]];
        }
        particles_ = resampled;
    }

    void apply_roughening()
    {
        if (roughening_scale_ <= Scalar{0}) return;

        Scalar np_factor = std::pow(static_cast<Scalar>(NP),
                                    Scalar{-1} / static_cast<Scalar>(NX));

        std::normal_distribution<Scalar> normal(Scalar{0}, Scalar{1});

        for (std::size_t d = 0; d < NX; ++d) {
            Scalar min_val = particles_[0](static_cast<int>(d));
            Scalar max_val = min_val;
            for (std::size_t i = 1; i < NP; ++i) {
                Scalar v = particles_[i](static_cast<int>(d));
                min_val = std::min(min_val, v);
                max_val = std::max(max_val, v);
            }
            Scalar spread = max_val - min_val;
            Scalar sigma = roughening_scale_ * spread * np_factor;

            // Fallback: if all particles collapsed, use process noise scale
            if (sigma <= Scalar{0}) {
                sigma = roughening_scale_ * std::sqrt(Q_(static_cast<int>(d),
                                                         static_cast<int>(d)));
            }

            if (sigma > Scalar{0}) {
                for (std::size_t i = 0; i < NP; ++i) {
                    particles_[i](static_cast<int>(d)) += sigma * normal(rng_);
                }
            }
        }
    }
};

// Factory function since NP cannot be deduced via CTAD
template<std::size_t NP, typename Dynamics, typename Measurement,
         typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY,
         typename Rng = std::mt19937_64>
auto make_particle_filter(Dynamics dynamics, Measurement measurement,
                          pf_config<Scalar, NX, NU, NY> config,
                          Rng rng = Rng{})
{
    return particle_filter<Scalar, NX, NU, NY, NP, Dynamics, Measurement,
                           systematic_resampling, Rng>(
        std::move(dynamics), std::move(measurement), std::move(config), std::move(rng));
}

// Static assert helpers
namespace detail {

struct pf_sa_dynamics {
    auto operator()(const Vector<double, 2>&,
                    const Vector<double, 1>&) const -> Vector<double, 2>
    {
        return Vector<double, 2>::Zero();
    }
};

struct pf_sa_measurement {
    auto operator()(const Vector<double, 2>&) const -> Vector<double, 1>
    {
        return Vector<double, 1>::Zero();
    }
};

using pf_test_type = particle_filter<double, 2, 1, 1, 10,
                                     pf_sa_dynamics, pf_sa_measurement>;

}

static_assert(ObserverPolicy<detail::pf_test_type>);
static_assert(!CovarianceObserver<detail::pf_test_type>);

}

#endif
