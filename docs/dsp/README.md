# Signal Processing

Digital signal processing primitives for filtering within control loops. These
are discrete-time filters designed for sample-by-sample processing in real-time
applications.

## Types

- [biquad](biquad.md) -- IIR second-order section with factory functions for low-pass, high-pass, band-pass, notch, and dirty derivative
- [fir](fir.md) -- Finite impulse response filter with compile-time tap count
- [discrete_filter](discrete-filter.md) -- Cascaded biquad filter for higher-order IIR responses

## When to use

Pick **biquad** for standard second-order filtering tasks: low-pass noise
rejection, notch filtering of known disturbance frequencies, or dirty derivative
computation.

Pick **fir** for linear-phase filtering or when you have specific impulse response
coefficients.

Pick **discrete_filter** when you need higher-order IIR responses built from
cascaded second-order sections.

## Theory

For the mathematical background, see [DSP Theory](../reference/dsp-theory.md).
