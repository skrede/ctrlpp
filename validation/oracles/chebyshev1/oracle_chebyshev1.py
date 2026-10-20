#!/usr/bin/env python3
"""Exact Chebyshev Type I low-pass design, the reference for the two Chebyshev
precision cases in tests/unit/dsp_hardening_test.cpp.

The script evaluates the design definition at 60 significant digits with the
Python standard library only (decimal), independent of ctrlpp's arithmetic:

  Omega   = tan(pi fc / fs)                         bilinear pre-warp
  eps     = sqrt(10^(r / 10) - 1)                   ripple factor
  v       = asinh(1 / eps) / N
  theta_k = (2k + 1) pi / (2N),  k = 0 .. N/2 - 1   section order of the library
  sigma_k = -sinh(v) sin(theta_k),  w_k = cosh(v) cos(theta_k)
  g1 = -2 sigma_k,  g0 = sigma_k^2 + w_k^2          prototype g0 / (p^2 + g1 p + g0)

and the bilinear map p = (z - 1) / (Omega (z + 1)) of each prototype section,
divided through by its leading denominator coefficient
den0 = 1 + g1 Omega + g0 Omega^2:

  b0 = b2 = g0 Omega^2 / den0,  b1 = 2 b0
  a1 = (2 g0 Omega^2 - 2) / den0,  a2 = (1 - g1 Omega + g0 Omega^2) / den0

Each unnormalized section has DC gain one; the cascade's DC target is
1 / sqrt(1 + eps^2), the gain at DC of an even-order design whose passband peak
is 0 dB. Before printing, the script asserts that every unnormalized section
has DC gain one, that the target and the normalized cascade's gain at fc, the
passband edge, both equal the ripple definition 10^(-r / 20), and that its pi
and tan(pi / 10) agree with closed forms.

Every value is printed as a C++ hex-float pair: hi is the double nearest the
exact value and lo the double nearest the remainder.

The inputs are exact binary values given as hex-float strings, so a float
specification such as 1e-4f passes in exactly:

  python3 oracle_chebyshev1.py                       both designs the test uses
  python3 oracle_chebyshev1.py ORDER FC FS RIPPLE    for example 4 0x1.9p+6 0x1.f4p+9 0x1p+0
"""
import sys
import struct
from decimal import Decimal, getcontext

getcontext().prec = 60
GUARD = Decimal(10) ** -55


def pi():
    """Series of the decimal module documentation recipe."""
    getcontext().prec += 2
    three = Decimal(3)
    lasts, t, s, n, na, d, da = 0, three, 3, 1, 0, 0, 24
    while s != lasts:
        lasts = s
        n, na = n + na, na + 8
        d, da = d + da, da + 32
        t = (t * n) / d
        s += t
    getcontext().prec -= 2
    return +s


def arctan_inverse(m):
    """arctan(1 / m) by its Taylor series, for the Machin cross-check of pi."""
    getcontext().prec += 4
    x = Decimal(1) / m
    x2 = x * x
    term, total, k = x, x, 1
    while True:
        term *= -x2
        step = term / (2 * k + 1)
        if abs(step) < Decimal(10) ** -(getcontext().prec + 2):
            break
        total += step
        k += 1
    getcontext().prec -= 4
    return +total


def sin_cos(x):
    """Taylor series of the decimal module documentation recipe, both at once."""
    getcontext().prec += 4
    x2 = x * x
    s, c = x, Decimal(1)
    st, ct, i = x, Decimal(1), 1
    while True:
        st *= -x2 / ((2 * i) * (2 * i + 1))
        ct *= -x2 / ((2 * i - 1) * (2 * i))
        if abs(st) + abs(ct) < Decimal(10) ** -(getcontext().prec + 2):
            break
        s += st
        c += ct
        i += 1
    getcontext().prec -= 4
    return +s, +c


def sinh_cosh(x):
    e = x.exp()
    return (e - 1 / e) / 2, (e + 1 / e) / 2


def asinh(y):
    return (y + (y * y + 1).sqrt()).ln()


def chebyshev1_sections(order, fc, fs, ripple):
    p = pi()
    s, c = sin_cos(p * fc / fs)
    omega = s / c
    eps = ((Decimal(10) ** (ripple / 10)) - 1).sqrt()
    v = asinh(1 / eps) / order
    sinh_v, cosh_v = sinh_cosh(v)
    sections = []
    for k in range(order // 2):
        sin_t, cos_t = sin_cos(p * (2 * k + 1) / (2 * order))
        sigma = -sinh_v * sin_t
        w = cosh_v * cos_t
        g1, g0 = -2 * sigma, sigma * sigma + w * w
        den0 = 1 + g1 * omega + g0 * omega * omega
        b0 = g0 * omega * omega / den0
        sections.append({
            "b0": b0,
            "a1": (2 * g0 * omega * omega - 2) / den0,
            "a2": (1 - g1 * omega + g0 * omega * omega) / den0,
        })
    return sections, 1 / (1 + eps * eps).sqrt()


def response(sections, z):
    """Cascade response at complex z = (re, im), each section with b1 = 2 b0, b2 = b0."""
    def mul(u, w):
        return (u[0] * w[0] - u[1] * w[1], u[0] * w[1] + u[1] * w[0])

    def div(u, w):
        den = w[0] * w[0] + w[1] * w[1]
        return ((u[0] * w[0] + u[1] * w[1]) / den, (u[1] * w[0] - u[0] * w[1]) / den)

    zi = div((Decimal(1), Decimal(0)), z)
    zi2 = mul(zi, zi)
    h = (Decimal(1), Decimal(0))
    for sec in sections:
        b0 = sec["b0"]
        num = (b0 + 2 * b0 * zi[0] + b0 * zi2[0], 2 * b0 * zi[1] + b0 * zi2[1])
        den = (1 + sec["a1"] * zi[0] + sec["a2"] * zi2[0], sec["a1"] * zi[1] + sec["a2"] * zi2[1])
        h = mul(h, div(num, den))
    return h


def self_check(fc, fs, ripple, sections, target):
    p = pi()
    machin = 16 * arctan_inverse(5) - 4 * arctan_inverse(239)
    assert abs(p - machin) < GUARD, "pi series disagrees with Machin's formula"
    s, c = sin_cos(p / 10)
    closed = (1 - 2 / Decimal(5).sqrt()).sqrt()
    assert abs(s / c - closed) < GUARD, "tan(pi / 10) disagrees with sqrt(1 - 2 / sqrt(5))"
    for sec in sections:
        dc = 4 * sec["b0"] / (1 + sec["a1"] + sec["a2"])
        assert abs(dc - 1) < GUARD, "an unnormalized section's DC gain is not one"
    normalized = [dict(sec) for sec in sections]
    normalized[0]["b0"] *= target
    s, c = sin_cos(2 * p * fc / fs)
    h = response(normalized, (c, s))
    gain = (h[0] * h[0] + h[1] * h[1]).sqrt()
    edge = Decimal(10) ** (-ripple / 20)
    assert abs(target - edge) < GUARD, "the DC target is not 10^(-r / 20)"
    assert abs(gain - edge) < GUARD, "the normalized cascade's gain at fc is not 10^(-r / 20)"


def hex_pair(x):
    hi = float(x)
    lo = float(x - Decimal(hi))
    return f"{{{hi.hex()}, {lo.hex()}}}"


def single_precision_hex(x):
    return struct.unpack("f", struct.pack("f", x))[0].hex()


def report(order, fc_hex, fs_hex, ripple_hex):
    fc, fs, ripple = (Decimal(float.fromhex(h)) for h in (fc_hex, fs_hex, ripple_hex))
    sections, target = chebyshev1_sections(order, fc, fs, ripple)
    self_check(fc, fs, ripple, sections, target)
    print(f"// order {order}, fc {fc_hex}, fs {fs_hex}, ripple {ripple_hex} dB")
    for k, sec in enumerate(sections):
        for name in ("a1", "a2", "b0"):
            print(f"section {k} {name}: {hex_pair(sec[name])}  // {sec[name]:.25e}")
    print(f"dc target: {hex_pair(target)}  // {target:.25e}")


if __name__ == "__main__":
    if len(sys.argv) == 5:
        report(int(sys.argv[1]), *sys.argv[2:])
    elif len(sys.argv) == 1:
        report(4, "0x1.9p+6", "0x1.f4p+9", "0x1p+0")
        report(4, "0x1.9p+6", "0x1.f4p+9", single_precision_hex(1e-4))
    else:
        sys.exit(__doc__)
