# Signal Processing

Digital signal processing primitives for filtering within control loops. These
are discrete-time filters designed for sample-by-sample processing in real-time
applications.

## Types

- [biquad](biquad.md)<br/> IIR second-order section with factory functions for a Butterworth low-pass, a notch, and a dirty derivative; `cascaded_biquad` with Butterworth and Chebyshev Type I designs; and the per-channel `vector_biquad` and `vector_cascaded_biquad`
- [fir](fir.md)<br/> Finite impulse response filter with compile-time tap count
- [discrete_filter](discrete-filter.md)<br/> C++20 concepts for scalar and vector digital filters

## When to use

Pick **biquad** for standard second-order filtering tasks: low-pass noise
rejection, notch filtering of known disturbance frequencies, or dirty derivative
computation. The library designs no high-pass or band-pass filter; such a section
is built from its own `biquad_coeffs`.

Pick **vector_biquad** (from vector_biquad.h) to filter each element of a
fixed-size vector through its own copy of one design.

Pick **fir** for linear-phase filtering or when you have specific impulse response
coefficients.

Use **cascaded_biquad** (from cascaded_biquad.h, which biquad.h includes) when you need higher-order IIR responses
built from cascaded second-order sections. The **discrete_filter** concept lets
you write generic code that accepts any scalar filter type.

The convenience header `<ctrlpp/dsp.h>` includes `fir.h`, `biquad.h` (and with it
`cascaded_biquad.h` and `dsp_types.h`), `vector_biquad.h` and `discrete_filter.h`.

## Theory

For the mathematical background, see [DSP Theory](../../background/dsp.md).
