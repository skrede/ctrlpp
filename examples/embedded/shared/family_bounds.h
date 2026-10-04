#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_FAMILY_BOUNDS_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_FAMILY_BOUNDS_H

#include <cstddef>

namespace ctrlpp {

/// @brief Roundings one predict-and-update step of a Kalman filter charges to
/// its covariance, at nx states and ny outputs, for a plant whose state matrix
/// entries each carry at most `plant_roundings`.
///
/// A product of a p-by-q and a q-by-r matrix is charged 2*p*q*r and a p-by-r
/// sum p*r. Predict forms A*P*A^T + Q. The gain forms S = C*P*C^T + R and C*P
/// and solves S against C*P by column-pivoting Householder QR, charged 4*m^3
/// with m = max(nx, ny) as the Riccati stages are. The Joseph update forms K*C,
/// I - K*C, the two products around P, K*R*K^T and their sum, and symmetrizes
/// with one add and one halving per entry. The plant's rounding reaches P
/// through both factors of A*P*A^T. The estimate, its innovation and the NIS do
/// not feed the covariance and are not charged; the Joseph form leaves the
/// covariance insensitive to the gain to first order, so the gain's charge is
/// slack.
///
/// The count charges every step without crediting the filter's forgetting. That
/// bounds the run when the first-order error map of the covariance is
/// non-expansive on average over it, which the generator evaluates on the
/// reference run before it prints a bound.
///
/// @cite verhaegen1986 -- Verhaegen & Van Dooren, "Numerical Aspects of Different Kalman Filter Implementations", IEEE Trans. Automatic Control 31(10), 1986 (first-order error propagation of the covariance recursion; the Joseph form's insensitivity to the gain)
constexpr double covariance_roundings_per_filter_step(std::size_t nx, std::size_t ny, double plant_roundings)
{
    const double n       = static_cast<double>(nx);
    const double p       = static_cast<double>(ny);
    const double m       = static_cast<double>(nx > ny ? nx : ny);
    const double predict = 4.0 * n * n * n + n * n;
    const double gain    = 4.0 * p * n * n + 2.0 * p * p * n + p * p + 4.0 * m * m * m;
    const double joseph  = 4.0 * n * n * n + 4.0 * n * n * p + 2.0 * n * p * p + 4.0 * n * n;
    return predict + gain + joseph + 2.0 * plant_roundings;
}

/// @brief Roundings one step of a transposed direct-form II low-pass biquad
/// charges at `biquad_scale`, the design's coefficient error included.
///
/// The step forms y = b0*x + w0 (2), w0' = b1*x - a1*y + w1 (4) and
/// w1' = b2*x - a2*y (3). A rounding of y perturbs w0 at the same step, so every
/// rounding lands on the state and reaches the output through the powers of
/// F = [-a1 1; -a2 0]; the worse state component collects y's and w0''s, 6.
///
/// The design forms w0 = 2*pi*fc/fs (3 roundings), its cosine and sine (each
/// one rounding for the library and, for w0 <= 1, the argument's three carried
/// through a slope of at most one), alpha = sin/sqrt(2) (2), 1 + alpha and its
/// reciprocal (2), and one or two more operations per coefficient. Measured at
/// its surrogate -- its formula in absolute values with |cos| <= 1 -- a1
/// carries 13 roundings, b0, b1 and b2 14, and a2 16. Those errors multiply
/// operands that sum to at most the scale, adding 16: 22 per step.
///
/// The count charges every step without crediting the filter's decay, which
/// bounds the run when the sum of the infinity norms of the powers of F over it
/// is at most its step count; the generator evaluates that before it prints a
/// bound.
///
/// @cite higham2002 -- Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., 2002, Ch. 3 (componentwise rounding of products and sums)
/// @cite bristowjohnson2005 -- Bristow-Johnson, "Cookbook Formulae for Audio EQ Biquad Filter Coefficients", 2005 (the low-pass formulas the design evaluates)
constexpr double biquad_roundings_per_step()
{
    return 6.0 + 16.0;
}

/// @brief Scale at which `biquad_roundings_per_step` charges, for a run whose
/// input and output never exceed `peak` in magnitude.
///
/// At their surrogates b0 and b2 are at most 1, b1 and a1 at most 2, and a2 is
/// 1, so the coefficients sum to c <= 7. The state follows from the output,
/// w0 = y - b0*x and w1 = w0' - b1*x + a1*y, so every sum the step forms is at
/// most (1 + 2*c) * peak.
constexpr double biquad_scale(double peak)
{
    return (1.0 + 2.0 * 7.0) * peak;
}

/// @brief Roundings one position sample of a rest-to-rest double-S profile
/// charges on the path that reaches its velocity limit.
///
/// Construction on that path is closed form: the displacement and its sign fold
/// (2), the two boundary velocities' sign folds (2), T_j = a/j (1),
/// T_a = T_j + v/a (2), T_v = h/v - T_a (2) and the total duration (2): 11. The
/// longest evaluation branch, the constant-deceleration segment, forms the
/// position in 27 and folds it back to the commanded frame in 3. The sample
/// time is an index times a rounded period: 2. A sample beside a segment
/// boundary may take the neighboring branch at one precision; the profile is
/// continuous there and no branch is longer, so the count covers either.
constexpr double double_s_roundings_per_sample()
{
    return 11.0 + 27.0 + 3.0 + 2.0;
}

/// @brief Scale at which every one of those roundings reaches the position.
///
/// Every position the evaluation forms is at most |q0| + |q1|, or one term of a
/// limit times a power of a local time no longer than the duration T. Every
/// time it forms is at most T, and the position's slope against any of them is
/// at most v + a*T + j*T^2. Each rounding therefore moves the sampled position
/// by at most epsilon times |q0| + |q1| + v*T + a*T^2 + j*T^3.
constexpr double double_s_scale(double q0, double q1, double v, double a, double j, double duration)
{
    const double t = duration;
    return (q0 < 0.0 ? -q0 : q0) + (q1 < 0.0 ? -q1 : q1) + v * t + a * t * t + j * t * t * t;
}

}

#endif
