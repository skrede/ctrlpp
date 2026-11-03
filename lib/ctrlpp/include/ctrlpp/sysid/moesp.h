#ifndef HPP_GUARD_CTRLPP_SYSID_MOESP_H
#define HPP_GUARD_CTRLPP_SYSID_MOESP_H

// @cite verhaegen2007 -- Verhaegen & Verdult, "Filtering and System Identification: A Least Squares Approach", 2007, Sec. 9.6 (PO-MOESP)

#include "ctrlpp/sysid/fit_metrics.h"
#include "ctrlpp/sysid/sysid_types.h"
#include "ctrlpp/sysid/sysid_result.h"

#include "ctrlpp/sysid/detail/moesp_steps.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/model/state_space.h"

#include <Eigen/Dense>

#include <cmath>
#include <cstddef>
#include <utility>
#include <algorithm>

namespace ctrlpp
{

namespace detail
{

template <typename Scalar>
struct moesp_subspace
{
    Scalar condition_number;
    Eigen::MatrixX<Scalar> gamma;
    Eigen::VectorX<Scalar> singular_values;
};

inline std::size_t resolve_block_rows(std::size_t block_rows, Eigen::Index samples)
{
    if(block_rows != 0)
        return block_rows;
    return static_cast<std::size_t>(std::min(samples / 4, Eigen::Index{30}));
}

template <typename Derived1, typename Derived2>
expected<void, sysid_error> check_moesp_record(const Eigen::MatrixBase<Derived1>& Y, const Eigen::MatrixBase<Derived2>& U, std::size_t i, std::size_t order)
{
    if(Y.cols() != U.cols())
        return unexpected(sysid_error::record_length_mismatch);
    if(Y.rows() != 1 || U.rows() != 1)
        return unexpected(sysid_error::record_not_single_row);
    auto const ny = static_cast<std::size_t>(Y.rows());
    auto const samples = static_cast<std::size_t>(Y.cols());
    if(i == 0 || (i - 1) * ny < order)
        return unexpected(sysid_error::block_rows_below_order);
    // N - 2i + 1 < i + order, rearranged so no unsigned subtraction wraps; i > N
    // already implies it and bounds 3i away from overflow.
    if(i > samples || samples + 1 < 3 * i + order)
        return unexpected(sysid_error::too_few_samples);
    if(!Y.allFinite() || !U.allFinite())
        return unexpected(sysid_error::non_finite_sample);
    return {};
}

// @cite verhaegen2007 -- Gamma = U_1 S_1^{1/2} from the SVD of the projection
template <std::size_t NX, typename Scalar>
expected<moesp_subspace<Scalar>, sysid_error> moesp_observability(const Eigen::MatrixX<Scalar>& O_i)
{
    static constexpr auto nx = static_cast<Eigen::Index>(NX);
    Eigen::BDCSVD<Eigen::MatrixX<Scalar>> svd(O_i, Eigen::ComputeThinU | Eigen::ComputeThinV);
    auto const& sv = svd.singularValues();
    auto const& U_svd = svd.matrixU();

    Scalar const cond = sv(0) / sv(nx - 1);
    if(!std::isfinite(cond))
        return unexpected(sysid_error::model_not_identified);

    Eigen::MatrixX<Scalar> Gamma(U_svd.rows(), nx);
    for(Eigen::Index c = 0; c < nx; ++c)
        Gamma.col(c) = U_svd.col(c) * std::sqrt(sv(c));
    if(!Gamma.allFinite())
        return unexpected(sysid_error::model_not_identified);

    return moesp_subspace<Scalar>{.condition_number = cond, .gamma = std::move(Gamma), .singular_values = sv};
}

// Filling the result here rather than in moesp_system keeps that function's frame small; the result reaches
// the caller through a single move.
template <typename Scalar, std::size_t NX>
moesp_result<Scalar, NX, 1, 1> assemble_moesp_result(
        const discrete_state_space<Scalar, NX, 1, 1>& sys, moesp_subspace<Scalar>& subspace, const fit_metrics<Scalar>& metrics, const Vector<Scalar, NX>& x0)
{
    moesp_result<Scalar, NX, 1, 1> result;
    result.system = sys;
    result.singular_values = std::move(subspace.singular_values);
    result.metrics = metrics;
    result.condition_number = subspace.condition_number;
    result.initial_state = x0;
    return result;
}

template <std::size_t NX, typename Scalar, typename Derived1, typename Derived2>
expected<moesp_result<Scalar, NX, 1, 1>, sysid_error> moesp_system(moesp_subspace<Scalar> subspace, const Eigen::MatrixBase<Derived1>& Y, const Eigen::MatrixBase<Derived2>& U)
{
    // The output row is the first block row of Gamma; its first column only
    // coincides with it for NX = 1.
    Matrix<Scalar, 1, NX> C = subspace.gamma.topRows(1);

    auto A = extract_system_A<Scalar, NX>(subspace.gamma, Y.rows());
    if(!A.allFinite())
        return unexpected(sysid_error::model_not_identified);

    auto const fit = recover_initial_state_and_BD<Scalar, NX>(A, C, Y, U);
    if(!fit.B.allFinite() || !fit.D.allFinite() || !fit.x0.allFinite())
        return unexpected(sysid_error::model_not_identified);

    discrete_state_space<Scalar, NX, 1, 1> sys{A, fit.B, C, fit.D};
    return assemble_moesp_result(sys, subspace, compute_moesp_metrics(sys, fit.x0, Y, U), fit.x0);
}

}

template <typename Derived1, typename Derived2>
expected<Eigen::VectorX<typename Derived1::Scalar>, sysid_error> moesp_singular_values(
        const Eigen::MatrixBase<Derived1>& Y, const Eigen::MatrixBase<Derived2>& U, std::size_t block_rows = 0)
{
    using Scalar = typename Derived1::Scalar;

    auto const i = detail::resolve_block_rows(block_rows, Y.cols());
    if(auto const gate = detail::check_moesp_record(Y, U, i, 1); !gate)
        return unexpected(gate.error());

    auto const projection = detail::moesp_projection(Y, U, static_cast<Eigen::Index>(i));
    Eigen::BDCSVD<Eigen::MatrixX<Scalar>> svd(projection.O_i, Eigen::ComputeThinU | Eigen::ComputeThinV);
    return Eigen::VectorX<Scalar>{svd.singularValues()};
}

template <std::size_t NX, typename Derived1, typename Derived2>
    requires ctrlpp_floating_scalar<typename Derived1::Scalar>
expected<moesp_result<typename Derived1::Scalar, NX, 1, 1>, sysid_error> moesp(const Eigen::MatrixBase<Derived1>& Y, const Eigen::MatrixBase<Derived2>& U, std::size_t block_rows = 0)
{
    static_assert(NX >= 1, "moesp requires NX >= 1");
    auto const i = detail::resolve_block_rows(block_rows, Y.cols());
    if(auto const gate = detail::check_moesp_record(Y, U, i, NX); !gate)
        return unexpected(gate.error());

    auto const projection = detail::moesp_projection(Y, U, static_cast<Eigen::Index>(i));
    auto subspace = detail::moesp_observability<NX>(projection.O_i);
    if(!subspace)
        return unexpected(subspace.error());
    return detail::moesp_system<NX>(std::move(*subspace), Y, U);
}

}

#endif
