#!/usr/bin/env python3
"""Bakeoff winner selection for the CARE method race.

Two-phase selection rule:

  Phase 1 (eligibility):  A candidate method `m` is eligible iff
    all(instr_median(m, nx) <= instr_median(ct_optcon, nx)
        for nx in [8, 12, 16, 20, 24, 30]).
    This is the strict CARE primary gate.

  Phase 2 (tie-break):    Among eligible candidates, compute
    total_instructions(m) = sum(instr_median(m, nx) for nx in TARGET_NX).
    The winner is the candidate with the minimum total_instructions.
    If two or more candidates are within 1% of the minimum total
    instruction count, secondary tie-break is
      total_wall_ns(m) = sum(wall_median(m, nx) for nx in TARGET_NX),
    argmin. This two-phase rule keeps the instruction gate primary
    while giving wall-clock the deciding vote when instructions are
    effectively tied.

If no candidate is eligible, the script prints `VERDICT: NO_WINNER`
to stderr and exits with code 2 so the caller can trigger the
conditional structured-Schur fallback.

The `--exclude <method>` flag removes a single candidate from the
METHOD_KEYS list before selection runs; used when a winner is
disqualified downstream (for example by a correctness test failure)
and a second-best candidate must be chosen mechanically without
hand-editing the script.

Usage:
  tools/bakeoff_winner.py --input <bench_care_methods.csv>
  tools/bakeoff_winner.py --input <bench_care_methods.csv> --exclude sign
"""
import argparse
import csv
import re
import sys

TARGET_NX   = [8, 12, 16, 20, 24, 30]
TIE_BAND    = 0.01  # 1% of minimum total instructions
NAME_RE     = re.compile(r'(?:ctrlpp::care\[(?P<method>\w+)\]|ct::optcon::CARE)\s+NX=(?P<nx>\d+)')
METHOD_KEYS = ['schur', 'sign', 'balanced']

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

def load(path):
    instr = {}
    wall  = {}
    for r in csv.DictReader(_open_nanobench_csv(path)):
        m = NAME_RE.search(r['name'])
        if not m:
            continue
        method = m.group('method') if m.group('method') else 'ct_optcon'
        nx     = int(m.group('nx'))
        instr[(method, nx)] = float(r['instructions'])
        wall [(method, nx)] = float(r['elapsed'])
    return instr, wall

def main():
    ap = argparse.ArgumentParser(
        description='Select the CARE bakeoff winner by total median instructions with a 1% tie-band falling back to total wall-clock. See module docstring for the two-phase rule.')
    ap.add_argument('--input', required=True, help='bench_care_methods.csv from the bakeoff run')
    ap.add_argument('--exclude', action='append', default=[],
                    help='exclude a METHOD_KEYS entry from selection; may be given multiple times')
    args = ap.parse_args()

    method_keys = [m for m in METHOD_KEYS if m not in set(args.exclude)]
    if not method_keys:
        print('ERROR: every candidate was excluded; refusing to select', file=sys.stderr)
        sys.exit(2)

    instr, wall = load(args.input)

    eligible = []
    for method in method_keys:
        try:
            passes = all(instr[(method, nx)] <= instr[('ct_optcon', nx)] for nx in TARGET_NX)
        except KeyError as e:
            print(f'SKIP method={method}: missing data {e}', file=sys.stderr)
            continue
        total_instr = sum(instr[(method, nx)] for nx in TARGET_NX)
        total_wall  = sum(wall [(method, nx)] for nx in TARGET_NX)
        print(f'method={method:10s}  passes_gate={passes}  total_instr={total_instr:.0f}  total_wall_ns={total_wall:.3e}')
        if passes:
            eligible.append((method, total_instr, total_wall))

    if not eligible:
        print('VERDICT: NO_WINNER (no candidate passes the strict <= gate at every NX in '
              f'{TARGET_NX}).', file=sys.stderr)
        sys.exit(2)

    min_instr = min(total for _, total, _ in eligible)
    tied = [(method, total, wall_ns) for method, total, wall_ns in eligible
            if (total - min_instr) / min_instr <= TIE_BAND]
    if len(tied) == 1:
        winner = tied[0][0]
    else:
        # Tie: break by total wall-clock argmin.
        winner = min(tied, key=lambda entry: entry[2])[0]
    print(f'VERDICT: WINNER={winner}')

if __name__ == '__main__':
    main()
