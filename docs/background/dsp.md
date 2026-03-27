# DSP Theory

Digital signal processing (DSP) in ctrlpp provides discrete-time filters for
signal conditioning, noise reduction, and derivative estimation. The library
implements IIR (biquad) and FIR filters, with cascading support for
higher-order responses.

## IIR vs FIR Filters

- **IIR** (Infinite Impulse Response): uses feedback (past outputs) in
  addition to past inputs. Achieves steep frequency rolloff with few
  coefficients but can be unstable if poorly designed.
- **FIR** (Finite Impulse Response): uses only past inputs (no feedback).
  Always stable, linear phase possible, but requires more coefficients for
  the same rolloff.

## Biquad Transfer Function

A biquad (second-order section) is the fundamental IIR building block. Its
z-domain transfer function is:

$$
H(z) = \frac{b_0 + b_1 z^{-1} + b_2 z^{-2}}{1 + a_1 z^{-1} + a_2 z^{-2}}
$$

Six coefficients ($b_0, b_1, b_2, a_1, a_2$, and an implicit $a_0 = 1$)
define any second-order transfer function. The corresponding difference
equation is:

$$
y[n] = b_0 \, x[n] + b_1 \, x[n-1] + b_2 \, x[n-2] - a_1 \, y[n-1] - a_2 \, y[n-2]
$$

ctrlpp uses the Direct Form II Transposed implementation for numerical
stability, requiring only two state variables per section.

### Common Biquad Types

- **Low-pass**: attenuates frequencies above a cutoff $f_c$
- **High-pass**: attenuates frequencies below a cutoff $f_c$
- **Band-pass**: passes a frequency band, attenuates outside
- **Notch**: rejects a specific frequency
- **Dirty derivative**: differentiator with built-in low-pass rolloff

## Bilinear Transform

The bilinear (Tustin) transform maps a continuous-time transfer function
$H(s)$ to discrete time $H(z)$ by the substitution:

$$
s = \frac{2}{T_s} \cdot \frac{z - 1}{z + 1}
$$

where $T_s$ is the sampling period. This mapping preserves stability (stable
continuous filters map to stable discrete filters) and maps the entire
$j\omega$ axis to the unit circle. Frequency warping must be accounted for
when specifying cutoff frequencies:

$$
\omega_a = \frac{2}{T_s} \tan\!\left(\frac{\omega_d T_s}{2}\right)
$$

where $\omega_d$ is the desired digital frequency and $\omega_a$ is the
pre-warped analogue frequency.

## FIR Convolution

An FIR filter of order $M$ computes the output as a weighted sum of the
current and past $M$ input samples:

$$
y[n] = \sum_{k=0}^{M} h_k \, x[n - k]
$$

where $h_0, h_1, \ldots, h_M$ are the filter coefficients (impulse response).
FIR filters are always stable and can achieve exact linear phase when the
coefficients are symmetric.

## Cascading

Higher-order filters are built by cascading biquad sections in series. The
overall transfer function is the product of the individual sections:

$$
H(z) = \prod_{i=1}^{L} H_i(z)
$$

A fourth-order Butterworth low-pass, for example, is two cascaded
second-order sections ($L = 2$). ctrlpp provides `cascaded_biquad` for this
pattern and a `make_butterworth` factory for common filter designs.

## Frequency Response

The frequency response of a discrete filter is evaluated on the unit circle
$z = e^{j \omega T_s}$:

$$
H(e^{j\omega T_s}) = \frac{b_0 + b_1 e^{-j\omega T_s} + b_2 e^{-2j\omega T_s}}{1 + a_1 e^{-j\omega T_s} + a_2 e^{-2j\omega T_s}}
$$

The magnitude $|H|$ gives the gain at each frequency; the phase
$\angle H$ gives the phase shift.

## References

- Oppenheim, A. V., Willsky, A. S., and Nawab, S. H. (1997). *Signals and
  Systems.* Prentice Hall, 2nd ed. [`oppenheim1997`]
  Standard textbook covering discrete-time systems, z-transform, frequency
  response, and filter design.

- Bristow-Johnson, R. (2005). "Cookbook Formulae for Audio EQ Biquad Filter
  Coefficients." [`bristowjohnson2005`]
  Widely referenced resource for biquad coefficient calculation (low-pass,
  high-pass, notch, peaking, shelving).

## Related API Pages

- [biquad](../api/dsp/biquad.md) -- IIR second-order section filter, cascaded
  biquad, Butterworth factory
- [fir](../api/dsp/fir.md) -- finite impulse response filter
- [discrete_filter](../api/dsp/discrete-filter.md) -- filter concept and
  cascading
