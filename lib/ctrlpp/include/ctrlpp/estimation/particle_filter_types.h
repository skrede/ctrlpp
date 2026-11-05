#ifndef HPP_GUARD_CTRLPP_ESTIMATION_PARTICLE_FILTER_TYPES_H
#define HPP_GUARD_CTRLPP_ESTIMATION_PARTICLE_FILTER_TYPES_H

#include "ctrlpp/types.h"

#include "ctrlpp/util/concepts.h"

#include <cstddef>
#include <concepts>

namespace ctrlpp
{

enum class extraction_method
{
    weighted_mean,
    map
};

enum class weight_representation
{
    log,
    linear
};

enum class pf_update_error
{
    non_finite_state,
    non_finite_measurement,
    non_finite_result
};

// Ordered by severity; the status latches, so it never returns to a lower rung.
enum class pf_health
{
    ok,
    weights_reset,
    resampling_rolled_back,
    non_finite_estimate
};

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct pf_config
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NY > 0, "Output dimension NY must be positive");
    Matrix<Scalar, NX, NX> Q{Matrix<Scalar, NX, NX>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Vector<Scalar, NX> x0{Vector<Scalar, NX>::Zero()};
    Matrix<Scalar, NX, NX> P0{Matrix<Scalar, NX, NX>::Identity()};
    Scalar ess_threshold{Scalar{-1}};
    Scalar roughening_scale{static_cast<Scalar>(0.2)};
    extraction_method extraction{extraction_method::weighted_mean};
    weight_representation weights{weight_representation::log};
};

namespace detail
{

template <ctrlpp_floating_scalar Scalar, std::size_t NY>
struct gaussian_likelihood
{
    Scalar operator()(const Vector<Scalar, NY>& z, const Vector<Scalar, NY>& z_pred, const Matrix<Scalar, NY, NY>& R_inv, Scalar log_det_2piR) const
    {
        Vector<Scalar, NY> innov = z - z_pred;
        Scalar mahal = (innov.transpose() * R_inv * innov)(0, 0);
        return Scalar{-0.5} * mahal - Scalar{0.5} * log_det_2piR;
    }
};

}

template <typename L, typename Scalar, std::size_t NY>
concept pf_likelihood_model = requires(const L& l, const Vector<Scalar, NY>& v, const Matrix<Scalar, NY, NY>& Rinv, Scalar s) {
    { l(v, v, Rinv, s) } -> std::convertible_to<Scalar>;
};

}

#endif
