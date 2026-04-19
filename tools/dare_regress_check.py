#!/usr/bin/env python3
"""D-14 DARE non-regression gate.

Instruction-count is the governor-independent primary metric: the symmetric
band abs(instr_delta) <= 1% flags both slowdowns (regression) and unexpected
speedups (algorithmic disturbance of a path we claim not to touch).

Wall-clock is noise-dominated on boost-capable CPUs inside the performance
governor; the gate enforces only the non-regression direction: a positive
wall_delta above 3% fails, and a negative wall_delta (measurement run is
faster than baseline) passes regardless of magnitude. This is strict
non-regression as documented in D-14, not symmetric.

Default target: NX in {8, 12, 16, 20, 24, 30} against the frozen baseline
2026-04-18_15-54_bench_dare_vs_ct_wide_sweep.csv.

Usage:  tools/dare_regress_check.py <baseline.csv> <current.csv> [--label ctrlpp::dare]
"""
import argparse
import csv
import sys

TARGET_NX = [8, 12, 16, 20, 24, 30]
INSTR_TOL = 0.01
WALL_TOL  = 0.03

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

def load_map(path, column):
    return {r['name']: float(r[column]) for r in csv.DictReader(_open_nanobench_csv(path))}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('baseline')
    ap.add_argument('current')
    ap.add_argument('--label', default='ctrlpp::dare')
    args = ap.parse_args()

    base_i = load_map(args.baseline, 'instructions')
    cur_i  = load_map(args.current,  'instructions')
    base_w = load_map(args.baseline, 'elapsed')
    cur_w  = load_map(args.current,  'elapsed')

    ok = True
    for nx in TARGET_NX:
        key = f'{args.label} NX={nx}'
        if key not in base_i or key not in cur_i:
            print(f'SKIP NX={nx}: missing {key}', file=sys.stderr)
            ok = False
            continue
        di = (cur_i[key] - base_i[key]) / base_i[key] if base_i[key] > 0 else float('inf')
        dw = (cur_w[key] - base_w[key]) / base_w[key] if base_w[key] > 0 else float('inf')
        instr_ok = abs(di) <= INSTR_TOL
        wall_ok  = dw <= WALL_TOL
        passed   = instr_ok and wall_ok
        print(f'{key:30s} instr_delta={di:+.2%} wall_delta={dw:+.2%} '
              f'{"OK" if passed else "FAIL"} (instr_band=+/-{INSTR_TOL:.0%}, wall_band<=+{WALL_TOL:.0%})')
        ok = ok and passed
    sys.exit(0 if ok else 1)

if __name__ == '__main__':
    main()
