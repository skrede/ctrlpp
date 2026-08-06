#!/usr/bin/env python3
"""Author a curated fuzz seed for the two algebraic Riccati randomized targets.

WHY THIS EXISTS. The curated corpus documentation requires every seed to be
written from named values by an encoder, because a hand-assembled seed is
unverifiable and undebuggable: nobody can say afterwards which field a byte
belonged to, and nobody can regenerate it when the decoder changes. This is that
encoder, committed so the corpus is regenerable by someone other than whoever
still holds a scratch script.

THE FIELD-ORDER CONTRACT, in the words the corpus provenance table carries.

    11 little-endian binary64, in decode order:
        A(0,0), A(0,1), A(1,0), A(1,1),   the state matrix by rows
        B(0), B(1),                       the input matrix
        Q(0,0), Q(0,1), Q(1,0), Q(1,1),   the raw weight factor by rows
        R                                 the raw input weight
    minimum 88 bytes

Both targets read exactly this order, so one encoder serves both and the order is
a single contract rather than two. Read it off each target's decoder before
trusting this comment; if a decoder gains, drops or reorders a field, every seed
is regenerated and both provenance tables are updated in the same change.

WHAT THE TARGET DOES TO THESE VALUES BEFORE IT USES THEM. Each of the first ten
fields is clamped to the bound the decoder states and then flushed to exact zero
below the floor the decoder derives from that bound, so a value below the floor
is NOT the value the target sees. The eleventh is clamped and then squared and
added to a derived floor, so the input weight the target solves with is never the
number written here. This tool reports the decoded configuration after every
write, which is what makes "check the emitted values against the target's own
clamps" a step the tool performs rather than an instruction a reader may skip.

Values are accepted by NAME and never by position, so a caller cannot silently
transpose two fields, and in hexadecimal float form as well as decimal, because
a recorded reproducing input is bit-specific and a decimal rendering of one does
not survive the round trip.

Usage:
    tools/fuzz_corpus_encode.py --out FILE \\
        --a00 V --a01 V --a10 V --a11 V --b0 V --b1 V \\
        --q00 V --q01 V --q10 V --q11 V --r V

    tools/fuzz_corpus_encode.py --decode FILE
"""
import sys
import struct
import argparse

# The decode order, and the only place it is written down in this file.
FIELD_ORDER = ['a00', 'a01', 'a10', 'a11', 'b0', 'b1', 'q00', 'q01', 'q10', 'q11', 'r']

# The decoder's own constants, mirrored here so the report below describes what
# the target will do rather than what this tool did. Both targets carry the same
# three: the clamp bound the entries are limited to, the zero floor one part in a
# hundred beneath it, and the input-weight floor derived from the bound squared.
CLAMP_BOUND = 2.0
ZERO_FLOOR = CLAMP_BOUND / 1e2
WEIGHT_FLOOR = (CLAMP_BOUND * CLAMP_BOUND) / 1e3

# The decoder reads eleven binary64 and refuses anything shorter.
MINIMUM_SIZE = 8 * len(FIELD_ORDER)


def parse_scalar(text):
    """A binary64 from decimal or from a hexadecimal float literal."""
    stripped = text.strip()
    body = stripped[1:] if stripped[:1] in '+-' else stripped
    if body[:2].lower() == '0x':
        return float.fromhex(stripped)
    return float(stripped)


def clamp_entry(value):
    """The clamp-then-flush the decoder applies to the ten matrix entries."""
    clamped = min(max(value, -CLAMP_BOUND), CLAMP_BOUND)
    return 0.0 if abs(clamped) < ZERO_FLOOR else clamped


def decoded_weight(value):
    """The input weight the target solves with, from the raw eleventh field."""
    clamped = min(max(value, -CLAMP_BOUND), CLAMP_BOUND)
    return clamped * clamped + WEIGHT_FLOOR


def encode(values):
    """The 88 bytes, little-endian, in decode order."""
    return struct.pack('<11d', *[values[name] for name in FIELD_ORDER])


def decode(payload):
    """The eleven raw fields a target reads back out of the first 88 bytes."""
    if len(payload) < MINIMUM_SIZE:
        raise ValueError(f'{len(payload)} bytes is below the decoder minimum of {MINIMUM_SIZE}')
    unpacked = struct.unpack('<11d', payload[:MINIMUM_SIZE])
    return dict(zip(FIELD_ORDER, unpacked))


def report(values, stream):
    """Every field as written, as the decoder will see it, and in hexadecimal."""
    stream.write(f'{"field":>6}  {"written":>24}  {"hexadecimal":>24}  {"target sees":>24}\n')
    for name in FIELD_ORDER:
        written = values[name]
        if name == 'r':
            seen = f'{decoded_weight(written):.17g} (squared + floor)'
        else:
            entry = clamp_entry(written)
            note = ''
            if entry == 0.0 and written != 0.0:
                note = ' (flushed to zero)'
            elif entry != written:
                note = ' (clamped to the bound)'
            seen = f'{entry:.17g}{note}'
        stream.write(f'{name:>6}  {written:>24.17g}  {written.hex():>24}  {seen}\n')


def main():
    parser = argparse.ArgumentParser(
        description='Author a curated seed for the algebraic Riccati randomized targets '
                    'from named values, or decode one that already exists.',
        epilog='Field order: ' + ', '.join(FIELD_ORDER) +
               '; 11 little-endian binary64; minimum %d bytes.' % MINIMUM_SIZE)
    parser.add_argument('--out', help='seed file to write')
    parser.add_argument('--decode', help='seed file to read back and report, writing nothing')
    for name in FIELD_ORDER:
        parser.add_argument('--' + name, help='decimal or hexadecimal float')
    args = parser.parse_args()

    if args.decode:
        with open(args.decode, 'rb') as handle:
            payload = handle.read()
        values = decode(payload)
        print(f'{args.decode}: {len(payload)} bytes')
        report(values, sys.stdout)
        return 0

    if not args.out:
        parser.error('one of --out or --decode is required')

    missing = [name for name in FIELD_ORDER if getattr(args, name) is None]
    if missing:
        parser.error('every field must be named: missing ' + ', '.join(missing))

    values = {name: parse_scalar(getattr(args, name)) for name in FIELD_ORDER}

    with open(args.out, 'wb') as handle:
        handle.write(encode(values))

    # The write is not trusted. The file is read back, decoded through the same
    # contract a target applies, and compared BITWISE against what was asked for;
    # a report alone would let a rounding in the parse pass unnoticed.
    with open(args.out, 'rb') as handle:
        payload = handle.read()
    readback = decode(payload)

    mismatched = [name for name in FIELD_ORDER
                  if struct.pack('<d', readback[name]) != struct.pack('<d', values[name])]
    if mismatched:
        print('READBACK MISMATCH on ' + ', '.join(mismatched), file=sys.stderr)
        return 1

    print(f'{args.out}: {len(payload)} bytes, readback bitwise identical')
    report(readback, sys.stdout)
    return 0


if __name__ == '__main__':
    sys.exit(main())
