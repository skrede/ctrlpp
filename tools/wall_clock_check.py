#!/usr/bin/env python3
"""D-13 secondary gate: median wall-clock must not regress against a baseline CSV.

Usage:  tools/wall_clock_check.py <baseline.csv> <current.csv> [--label-substring ctrlpp::] [--tolerance 0.0]

By default the gate is strict (delta <= 0) per D-13. A --tolerance flag
allows a small measurement-noise band if the operator chooses (CONTEXT.md
D-13 is silent on the exact noise band; default strict).
"""
import argparse
import csv
import sys

def _open_nanobench_csv(path):
    """Open a nanobench CSV, skipping any leading blank lines before the header."""
    fh = open(path)
    while True:
        pos = fh.tell()
        line = fh.readline()
        if not line or line.strip():
            fh.seek(pos)
            break
    return fh

def load_map(path, column='elapsed'):
    return {r['name']: float(r[column]) for r in csv.DictReader(_open_nanobench_csv(path))}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('baseline')
    ap.add_argument('current')
    ap.add_argument('--label-substring', default='ctrlpp::')
    ap.add_argument('--tolerance', type=float, default=0.0,
                    help='maximum fractional increase permitted (default 0.0 = strict non-regression)')
    args = ap.parse_args()

    base = load_map(args.baseline)
    cur  = load_map(args.current)
    ok = True
    for name, t_cur in sorted(cur.items()):
        if args.label_substring and args.label_substring not in name:
            continue
        if name not in base:
            continue
        t_base = base[name]
        delta = (t_cur - t_base) / t_base if t_base > 0 else float('inf')
        passed = delta <= args.tolerance
        print(f'{name:50s} base={t_base:.3e} cur={t_cur:.3e} delta={delta:+.2%} {"OK" if passed else "FAIL"}')
        ok = ok and passed
    sys.exit(0 if ok else 1)

if __name__ == '__main__':
    main()
