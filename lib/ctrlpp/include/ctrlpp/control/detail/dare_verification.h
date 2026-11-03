#ifndef HPP_GUARD_CTRLPP_CONTROL_DETAIL_DARE_VERIFICATION_H
#define HPP_GUARD_CTRLPP_CONTROL_DETAIL_DARE_VERIFICATION_H

#include "ctrlpp/control/dare_types.h"

#include "ctrlpp/detail/psd_acceptance.h"
#include "ctrlpp/detail/riccati_solution.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <limits>
#include <complex>
#include <cstddef>
#include <algorithm>

namespace ctrlpp::detail
{

// unresolved and gain_unavailable are both absences of evidence, neither a certificate nor a refutation. They are kept
// apart because only gain_unavailable has a known cause: R + B'PB left the range at the scale being checked.
enum class dare_verification
{
    verified,
    refuted,
    unresolved,
    gain_unavailable,
};

template <typename Scalar, std::size_t NX, std::size_t NU>
dare_verification compute_dare_gain(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& P, Eigen::Matrix<Scalar, int(NU), int(NX)>& gain)
{
    const auto BtP = (B.transpose() * P).eval();
    const auto gain_operand = (R + BtP * B).eval();
    if(!BtP.allFinite() || !gain_operand.allFinite())
        return dare_verification::gain_unavailable;

    auto gain_qr = gain_operand.colPivHouseholderQr();
    gain_qr.setThreshold(Scalar{static_cast<int>(NU)} * std::numeric_limits<Scalar>::epsilon());
    if(!gain_qr.isInvertible())
        return dare_verification::refuted;

    // A non-finite gain from a finite, invertible operand says the gain does not exist at this scale, as an operand
    // leaving the range does; neither says whether P solves the equation.
    gain = gain_qr.solve(BtP * A).eval();
    return gain.allFinite() ? dare_verification::verified : dare_verification::gain_unavailable;
}

template <typename Scalar, std::size_t NX, std::size_t NU>
bool dare_residual_into(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B, const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& P, const Eigen::Matrix<Scalar, int(NU), int(NX)>& gain, Eigen::Matrix<Scalar, int(NX), int(NX)>& residual)
{
    using MatNxN = Eigen::Matrix<Scalar, int(NX), int(NX)>;

    const MatNxN AtPA = (A.transpose() * P * A).eval();
    const MatNxN AtPBK = (A.transpose() * P * B * gain).eval();
    residual = (AtPA - P - AtPBK + Q).eval();
    return AtPA.allFinite() && AtPBK.allFinite() && residual.allFinite();
}

template <typename Scalar, std::size_t NX>
dare_verification dare_closed_loop_spectrum(const Eigen::Matrix<Scalar, int(NX), int(NX)>& closed_loop)
{
    constexpr int n = static_cast<int>(NX);

    Eigen::EigenSolver<Eigen::Matrix<Scalar, n, n>> eigensolver(closed_loop, false);
    if(eigensolver.info() != Eigen::Success || !eigensolver.eigenvalues().allFinite())
        return dare_verification::unresolved;

    const Scalar closed_loop_scale = std::max(Scalar{1}, closed_loop.cwiseAbs().maxCoeff());
    // The standard n eps ||A - BK|| backward-error margin of each closed-loop eigenvalue.
    const Scalar unit_margin = Scalar{static_cast<int>(NX)} * std::numeric_limits<Scalar>::epsilon() * closed_loop_scale;
    for(Eigen::Index index = 0; index < n; ++index)
    {
        if(!(std::abs(eigensolver.eigenvalues()(index)) < Scalar{1} - unit_margin))
            return dare_verification::refuted;
    }
    return dare_verification::verified;
}

template <typename Scalar, std::size_t NX, std::size_t NU>
dare_verification verify_dare_closed_loop(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NX), int(NX)>& P, const Eigen::Matrix<Scalar, int(NU), int(NX)>& gain)
{
    using MatNxN = Eigen::Matrix<Scalar, int(NX), int(NX)>;

    MatNxN residual;
    if(!dare_residual_into<Scalar, NX, NU>(A, B, Q, P, gain, residual))
        return dare_verification::unresolved;

    const MatNxN closed_loop = (A - B * gain).eval();
    if(!closed_loop.allFinite())
        return dare_verification::unresolved;

    // The acceptance quantity is the answer's estimated forward error, not the residual: no threshold on the residual
    // separates the answers that keep half the significand from those that do not.
    switch(riccati_forward_error_verdict<Scalar, int(NX)>(closed_loop, residual, P))
    {
        case riccati_accuracy::within:
            break;
        case riccati_accuracy::exceeded:
            return dare_verification::refuted;
        case riccati_accuracy::unresolved:
            return dare_verification::unresolved;
    }
    return dare_closed_loop_spectrum<Scalar, NX>(closed_loop);
}

// verified_gain is written only on a verified verdict.
template <typename Scalar, std::size_t NX, std::size_t NU>
dare_verification verify_dare_solution(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R, const Eigen::Matrix<Scalar, int(NX), int(NX)>& P,
        Eigen::Matrix<Scalar, int(NU), int(NX)>& verified_gain)
{
    if(!P.allFinite())
        return dare_verification::unresolved;

    Eigen::Matrix<Scalar, int(NU), int(NX)> gain;
    const dare_verification gain_state = compute_dare_gain<Scalar, NX, NU>(A, B, R, P, gain);
    if(gain_state != dare_verification::verified)
        return gain_state;

    // At the caller's scale P is an exact power-of-two multiple of the matrix the extraction checked, so this repeats
    // that verdict except where an entry left the normal range in the rescale.
    if(!psd_within_rounding<Scalar, int(NX)>(P))
        return dare_verification::refuted;

    const dare_verification verdict = verify_dare_closed_loop<Scalar, NX, NU>(A, B, Q, P, gain);
    if(verdict == dare_verification::verified)
        verified_gain = gain;
    return verdict;
}

template <typename Scalar, std::size_t NX, std::size_t NU>
dare_verification verify_dare_solution(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
        const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R, const Eigen::Matrix<Scalar, int(NX), int(NX)>& P)
{
    Eigen::Matrix<Scalar, int(NU), int(NX)> gain;
    return verify_dare_solution<Scalar, NX, NU>(A, B, Q, R, P, gain);
}

// Two posings of the same problem share one gain, so the gains formed at both scales are held to each other.
template <typename Scalar, std::size_t NX, std::size_t NU>
bool dare_gains_agree(const Eigen::Matrix<Scalar, int(NU), int(NX)>& scaled_gain, const Eigen::Matrix<Scalar, int(NU), int(NX)>& returned_gain)
{
    const resolved_magnitude<Scalar> gain_scale = largest_magnitude<Scalar>({resolve_magnitude(scaled_gain), resolve_magnitude(returned_gain)});
    const resolved_magnitude<Scalar> gain_margin = scaled_magnitude(gain_scale, dare_gain_agreement_margin<Scalar, NX, NU>());
    return magnitude_within(resolve_magnitude((scaled_gain - returned_gain).eval()), gain_margin);
}

// The caller-scale claim is carried from the equilibrated one, since the equation is homogeneous of degree one in
// (P, Q, R); what homogeneity cannot supply, finiteness and definiteness after the rescale, is checked. The direct check
// then declines on every verdict but verified and gain_unavailable, the one absence of evidence with a known cause.
template <typename Scalar, std::size_t NX, std::size_t NU>
bool carry_dare_solution(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B, const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
        const Eigen::Matrix<Scalar, int(NU), int(NU)>& R, Scalar weight_scale, dare_result<Scalar, NX, NU>& result)
{
    result.P *= weight_scale;
    if(!result.P.allFinite())
        return false;

    if(!psd_within_rounding<Scalar, int(NX)>(result.P))
        return false;

    Eigen::Matrix<Scalar, int(NU), int(NX)> returned_gain;
    const dare_verification direct = verify_dare_solution<Scalar, NX, NU>(A, B, Q, R, result.P, returned_gain);
    if(direct != dare_verification::verified && direct != dare_verification::gain_unavailable)
        return false;

    return direct != dare_verification::verified || dare_gains_agree<Scalar, NX, NU>(result.K, returned_gain);
}

}

#endif
