#!/usr/bin/env python3
"""D-12 primary gate: median instructions <= ct::optcon at NX in {8, 12, 16, 20, 24, 30}.

Usage:  tools/care_primary_check.py <bench_care_vs_ct.csv>  [--bench-label ctrlpp::care]
        tools/care_primary_check.py <bench_lqr_continuous_vs_ct.csv> --bench-label ctrlpp::lqr_gain_continuous --ct-label ct::optcon::LQR
"""
import argparse
import csv
import sys

TARGET_NX = [8, 12, 16, 20, 24, 30]

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

def load_csv(path):
    return {r['name']: r for r in csv.DictReader(_open_nanobench_csv(path))}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('csv', help='nanobench CSV output')
    ap.add_argument('--bench-label', default='ctrlpp::care')
    ap.add_argument('--ct-label', default='ct::optcon::CARE')
    args = ap.parse_args()

    rows = load_csv(args.csv)
    ok = True
    for nx in TARGET_NX:
        c_key = f'{args.bench_label} NX={nx}'
        r_key = f'{args.ct_label} NX={nx}'
        if c_key not in rows or r_key not in rows:
            print(f'SKIP NX={nx}: missing row (candidate={c_key in rows}, ct={r_key in rows})', file=sys.stderr)
            ok = False
            continue
        c = float(rows[c_key]['instructions'])
        r = float(rows[r_key]['instructions'])
        passed = c <= r
        delta = (c - r) / r if r > 0 else float('inf')
        print(f'NX={nx:2d}  ctrlpp={c:>14.0f}  ct={r:>14.0f}  delta={delta:+.2%}  {"OK" if passed else "FAIL"}')
        ok = ok and passed
    sys.exit(0 if ok else 1)

if __name__ == '__main__':
    main()
