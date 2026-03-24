# DSP Theory

Digital signal processing (DSP) in ctrlpp provides discrete-time filters for
signal conditioning, noise reduction, and derivative estimation. The library
implements IIR (biquad) and FIR filters, with cascading support for
higher-order responses.

## Key Concepts

### IIR vs FIR Filters

- **IIR** (Infinite Impulse Response): uses feedback (past outputs) in addition
  to past inputs. Achieves steep frequency rolloff with few coefficients but
  can be unstable if poorly designed.
- **FIR** (Finite Impulse Response): uses only past inputs (no feedback).
  Always stable, linear phase possible, but requires more coefficients for
  the same rolloff.

### Biquad Sections

A biquad (second-order section) is the fundamental IIR building block:

```
y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2]
```

Six coefficients (b0, b1, b2, a1, a2, and an implicit a0=1) define any
second-order transfer function. ctrlpp uses the direct form II transposed
implementation for numerical stability.

Common biquad types:

- **Low-pass**: attenuates frequencies above a cutoff
- **High-pass**: attenuates frequencies below a cutoff
- **Band-pass**: passes a frequency band, attenuates outside
- **Notch**: rejects a specific frequency
- **Dirty derivative**: differentiator with built-in low-pass rolloff

### Cascading

Higher-order filters are built by cascading biquad sections in series. A
fourth-order Butterworth low-pass, for example, is two cascaded second-order
sections. ctrlpp provides `cascaded_biquad` for this pattern and a
`make_butterworth` factory for common filter designs.

### Frequency Response

The frequency response H(z) of a discrete filter is evaluated on the unit
circle z = exp(j*omega*T):

```
H(exp(j*omega*T)) = (b0 + b1*z^{-1} + b2*z^{-2}) /
                     (1  + a1*z^{-1} + a2*z^{-2})
```

The magnitude |H| gives the gain at each frequency; the phase angle(H) gives
the phase shift.

### Bilinear Transform

The bilinear (Tustin) transform maps a continuous-time transfer function to
discrete time:

```
s = (2/T) * (z - 1) / (z + 1)
```

This preserves stability (stable continuous filters map to stable discrete
filters) and maps the entire jw axis to the unit circle. Frequency warping
must be accounted for when specifying cutoff frequencies.

## References

- **Oppenheim, A. V., Willsky, A. S., and Nawab, S. H.** *Signals and
  Systems.* Prentice Hall, 2nd ed., 1997. ISBN 978-0-13-814757-0.
  Standard textbook covering discrete-time systems, z-transform, frequency
  response, and filter design.

- **Bristow-Johnson, R.** "Cookbook Formulae for Audio EQ Biquad Filter
  Coefficients." Audio EQ Cookbook, 2005.
  Widely referenced resource for biquad coefficient calculation (low-pass,
  high-pass, notch, peaking, shelving).

## Related API Pages

- [biquad](../dsp/biquad.md) -- IIR second-order section filter, cascaded
  biquad, Butterworth factory
- [fir](../dsp/fir.md) -- finite impulse response filter
- [discrete_filter](../dsp/discrete-filter.md) -- filter concept and
  cascading
