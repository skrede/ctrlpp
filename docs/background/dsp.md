# Digital Signal Processing

Digital signal processing (DSP) encompasses the mathematical tools and
algorithms for manipulating discrete-time signals. In control systems, DSP
techniques are used for noise filtering, derivative estimation, and signal
conditioning of sensor measurements before they enter the control loop
[1, Ch. 1, pp. 1--20].

This page covers the fundamental building blocks: the Z-transform, IIR and
FIR filters, the biquad structure, and the bilinear transform for converting
continuous-time filter designs to discrete time.

## Discrete-Time Signals

A discrete-time signal $x[n]$ is a sequence of values indexed by integer $n$.
The fundamental difference from continuous-time signals is that information
exists only at the sample instants $n T_s$, where $T_s$ is the sampling
period and $f_s = 1/T_s$ is the sampling frequency
[1, Ch. 2, pp. 69--78].

The Nyquist-Shannon sampling theorem states that a band-limited signal can be
perfectly reconstructed from its samples if $f_s > 2 f_{\max}$, where
$f_{\max}$ is the highest frequency component [1, Sec. 7.2, pp. 514--520].

## Z-Transform

The Z-transform maps a discrete-time signal to a complex-valued function of
$z$, analogous to the Laplace transform for continuous-time signals
[1, Ch. 10, pp. 732--770]:

$$
X(z) = \sum_{n=0}^{\infty} x[n] \, z^{-n}
$$

Key properties [1, Sec. 10.3, pp. 746--752]:

- **Linearity**: $a \cdot x_1[n] + b \cdot x_2[n] \longleftrightarrow a X_1(z) + b X_2(z)$
- **Time delay**: $x[n-k] \longleftrightarrow z^{-k} X(z)$
- **Convolution**: $x_1[n] * x_2[n] \longleftrightarrow X_1(z) \cdot X_2(z)$

The transfer function of a discrete linear system relates the Z-transforms of
output and input:

$$
H(z) = \frac{Y(z)}{X(z)}
$$

Stability requires all poles of $H(z)$ to lie strictly inside the unit circle
$|z| < 1$.

## IIR Filters

Infinite Impulse Response (IIR) filters use both feedforward (past inputs)
and feedback (past outputs) terms. A general IIR filter of order $N$ has
transfer function [1, Ch. 10, pp. 770--780]:

$$
H(z) = \frac{\sum_{k=0}^{M} b_k z^{-k}}{1 + \sum_{k=1}^{N} a_k z^{-k}}
$$

with corresponding difference equation:

$$
y[n] = \sum_{k=0}^{M} b_k \, x[n-k] - \sum_{k=1}^{N} a_k \, y[n-k]
$$

IIR filters achieve steep frequency rolloff with few coefficients but can be
unstable if the poles move outside the unit circle due to coefficient
quantisation or poor design.

## Biquad (Second-Order Section)

The biquad is the fundamental IIR building block. Its transfer function is
[2]:

$$
H(z) = \frac{b_0 + b_1 z^{-1} + b_2 z^{-2}}{1 + a_1 z^{-1} + a_2 z^{-2}}
$$

Six coefficients ($b_0, b_1, b_2, a_1, a_2$, with implicit $a_0 = 1$) define
any second-order filter. The difference equation is:

$$
y[n] = b_0 \, x[n] + b_1 \, x[n-1] + b_2 \, x[n-2] - a_1 \, y[n-1] - a_2 \, y[n-2]
$$

### Direct Form II Transposed

The Direct Form II Transposed (DF2T) implementation uses only two state
variables per section, providing good numerical properties
[1, Sec. 10.6, pp. 780--785]:

$$
y[n] = b_0 \, x[n] + s_1[n]
$$

$$
s_1[n+1] = b_1 \, x[n] - a_1 \, y[n] + s_2[n]
$$

$$
s_2[n+1] = b_2 \, x[n] - a_2 \, y[n]
$$

### Common Biquad Types

Standard biquad coefficient formulas [2]:

- **Low-pass**: attenuates frequencies above cutoff $f_c$ with quality
  factor $Q$ controlling the resonance peak
- **High-pass**: attenuates frequencies below cutoff $f_c$
- **Band-pass**: passes a frequency band centred at $f_0$ with bandwidth
  determined by $Q$
- **Notch (band-reject)**: rejects a specific frequency $f_0$
- **Dirty derivative**: differentiator with built-in low-pass rolloff,
  providing derivative estimation with noise attenuation

## Bilinear Transform

The bilinear (Tustin) transform maps a continuous-time transfer function
$H(s)$ to discrete time $H(z)$ via the substitution
[1, Sec. 10.7, pp. 790--798]:

$$
s = \frac{2}{T_s} \cdot \frac{z - 1}{z + 1}
$$

This mapping preserves stability (stable continuous poles map to stable
discrete poles) and maps the entire imaginary axis to the unit circle.

### Frequency Warping

The bilinear transform introduces a nonlinear frequency mapping between
analogue frequency $\omega_a$ and digital frequency $\omega_d$:

$$
\omega_a = \frac{2}{T_s} \tan\!\left(\frac{\omega_d T_s}{2}\right)
$$

Pre-warping compensates for this distortion: given a desired digital cutoff
frequency $\omega_d$, compute the pre-warped analogue frequency $\omega_a$
and design the analogue prototype at $\omega_a$.

## FIR Filters

Finite Impulse Response (FIR) filters use only feedforward terms (no
feedback). An FIR filter of order $M$ computes
[1, Ch. 10, pp. 768--770]:

$$
y[n] = \sum_{k=0}^{M} h_k \, x[n-k]
$$

where $h_0, h_1, \ldots, h_M$ are the filter coefficients (impulse response).

FIR advantages:
- **Always stable**: no poles, no risk of instability
- **Linear phase**: achievable with symmetric coefficients ($h_k = h_{M-k}$)
- **Simple design**: window method, Parks-McClellan algorithm

FIR disadvantages:
- **More coefficients needed**: for the same frequency selectivity as an
  IIR filter, typically 5--10 times more coefficients

## Cascading

Higher-order filters are built by cascading biquad sections in series
[1, Sec. 10.6, p. 785]:

$$
H(z) = \prod_{i=1}^{L} H_i(z)
$$

A fourth-order Butterworth low-pass filter, for example, is two cascaded
second-order sections ($L = 2$). Cascading is numerically more robust than
implementing a single high-order transfer function, because each second-order
section has well-conditioned coefficients.

## Frequency Response

The frequency response of a discrete filter is obtained by evaluating on the
unit circle $z = e^{j\omega T_s}$ [1, Sec. 10.5, pp. 776--780]:

$$
H(e^{j\omega T_s}) = \frac{b_0 + b_1 e^{-j\omega T_s} + b_2 e^{-2j\omega T_s}}{1 + a_1 e^{-j\omega T_s} + a_2 e^{-2j\omega T_s}}
$$

The magnitude $|H(e^{j\omega T_s})|$ gives the gain at each frequency; the
phase $\angle H(e^{j\omega T_s})$ gives the phase shift. The frequency
response is periodic with period $f_s$ and symmetric about $f_s/2$.

## References

[1] A. V. Oppenheim, A. S. Willsky, and S. H. Nawab, "Signals and Systems,"
2nd ed., Prentice Hall, 1997.

[2] R. Bristow-Johnson, "Cookbook Formulae for Audio EQ Biquad Filter
Coefficients," 2005.

[3] N. S. Nise, "Control Systems Engineering," 7th ed., Wiley, 2015.
