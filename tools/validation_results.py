#!/usr/bin/env python3
"""Generate the results table of the validation page, and gate the run's not-run cases.

The input is the test runner's native results file and not its JUnit serialization:
the JUnit form reports a failing unbuilt test inside a skipped element with a zero
failure count, which would publish a waiver nobody wrote.

The gate refuses every not-run case that is not a disabled test carrying a waiver
file. The generator emits one row per registered case, in the order the results file
lists them; nothing in the region varies between two runs of the same inputs, so the
check mode can compare byte for byte.

Usage: tools/validation_results.py (--gate | --check PAGE | --write PAGE | --self-test)
                                   [--build-dir DIR] [--validation-dir DIR] [--results FILE]
"""
import argparse
import difflib
import os
import sys
import tempfile
import xml.etree.ElementTree as ET

BEGIN_MARKER = '<!-- results:begin: generated from a harness run; edit the harness, not this region -->'
END_MARKER = '<!-- results:end -->'
HEADER = ['| Component | Octave function | Worst-case digits | Max abs error | Verdict |',
          '|-----------|----------------|:-----------------:|:-------------:|:-------:|']
DISABLED = 'Disabled'
IGNORED_SIGNALS = ('time', 'step', 'sample')
NO_NUMBER = 'n/a'

# The two label columns the numbers cannot supply. Hand-maintained, and policed
# against the run in both directions so a new case cannot silently lose its row.
CASE_LABELS = {
    'dare_solution': ('`dare`', '`dare()`'),
    'care_solution': ('`care`', '`care()`'),
    'lqr_step_response': ('`lqr` (infinite horizon)', '`dlqr()`'),
    'lqr_finite_horizon': ('`lqr` (finite horizon)', 'backward Riccati'),
    'lqi_tracking': ('`lqi`', '`dlqr()` augmented'),
    'pid_linear_step': ('`pid` (linear PI)', '`lsim()`'),
    'pole_placement': ('`place`', '`place()`'),
    'c2d_zoh': ('`discretize` (ZOH)', '`c2d()`'),
    'analysis_poles': ('`analysis` (poles)', '`pole()`, `ctrb()`, `obsv()`'),
    'tf_ss_conversion': ('`tf2ss` / `ss2tf`', '`tf2ss()`, `ss2tf()`'),
    'kalman_tracking': ('`kalman_filter`', 'time-varying KF'),
    'luenberger_observer': ('`luenberger_observer`', '`place()` on dual'),
    'butterworth_filter': ('`butterworth` (4th order)', '`butter()`, `filter()`'),
    'fir_filter': ('`fir`', '`filter()`'),
    'cubic_spline_natural': ('`cubic_spline` (natural)', '`csape()`, `ppder()`'),
    'so3_exp_log': ('`so3::exp` / `so3::log`', '`rot2q()`, `q2rot()`'),
    'batch_arx_ident': ('`batch_arx`', '`arx()`'),
    'moesp_ident': ('`moesp` (cross-algorithm)', '`n4sid()`'),
}

RESULTS_TEMPLATE = ('<Site><Testing><TestList><Test>./{name}</Test></TestList>'
                    '<Test Status="{status}"><Name>{name}</Name><Results><NamedMeasurement'
                    ' name="Completion Status"><Value>{completion}</Value>'
                    '</NamedMeasurement></Results></Test></Testing></Site>')

# One synthetic input per state the gate can reach, and the exit status it owes
# each. The last state is what isolates the completion-status guard: every other
# not-run state also lacks a waiver file, so without it the two guards are only
# ever proven together and either one could be removed unnoticed.
GATE_STATES = (
    ('ran and agreed', 'passed', 'Completed', True, 0),
    ('ran and failed', 'failed', 'Failed', False, 1),
    ('waived, waiver file present', 'notrun', DISABLED, True, 0),
    ('waived, waiver file absent', 'notrun', DISABLED, False, 1),
    ('registered but never built', 'notrun', 'Required Files Missing', False, 1),
    ('test command absent', 'notrun', 'Unable to find executable', False, 1),
    ('skipped at runtime', 'notrun', 'SKIP_REGULAR_EXPRESSION_MATCHED', False, 1),
    ('never built, waiver present', 'notrun', 'Required Files Missing', True, 1),
)

def find_results(build_dir):
    """Locate the results file of the most recent run under a build directory."""
    testing = os.path.join(build_dir, 'Testing')
    if not os.path.isdir(testing):
        sys.exit(f'{testing}: no results directory; run the harness first')
    return os.path.join(testing, open(os.path.join(testing, 'TAG')).readline().strip(), 'Test.xml')

def load_entries(path):
    """Return (name, status, completion status) per case, in registration order."""
    testing = ET.parse(path).getroot().find('Testing')
    outcomes = {}
    for test in testing.findall('Test'):
        completion = [m.findtext('Value') for m in test.iter('NamedMeasurement')
                      if m.get('name') == 'Completion Status']
        outcomes[test.findtext('Name')] = (test.get('Status'), completion[0].strip())
    names = [entry.text.rsplit('/', 1)[-1] for entry in testing.find('TestList')]
    return [(name, *outcomes[name]) for name in names if name in outcomes]

def waiver_path(cases_dir, name):
    return os.path.join(cases_dir, name, 'waiver.cfg')

def waiver_reason(path):
    """Return a waiver's reason verbatim, with the format's optional quotes stripped."""
    for line in open(path):
        if line.startswith('WAIVER_REASON='):
            return line.split('=', 1)[1].strip().strip('"')
    return ''

def run_gate(entries, cases_dir):
    """The run's accounting: every not-run case must be a waived one, and none may disagree."""
    failures = []
    for name, status, completion in entries:
        if status == 'failed':
            failures.append(f'{name}: the run reports it failed')
        elif status != 'notrun':
            continue
        elif completion != DISABLED:
            failures.append(f'{name}: did not run and is not waived, '
                            f'and its completion status was "{completion}"')
        elif not os.path.exists(waiver_path(cases_dir, name)):
            failures.append(f'{name}: is a disabled test carrying no waiver file at '
                            f'{waiver_path(cases_dir, name)}')
    for message in failures:
        print(message, file=sys.stderr)
    return 1 if failures else 0

def read_report(path):
    """Return a report's verdict and the digits and error of its least-agreeing signal."""
    verdict, worst = '', None
    for line in open(path):
        if '**Verdict**' in line:
            verdict = 'PASS' if '**PASS**' in line else 'FAIL'
            continue
        cells = [cell.strip() for cell in line.strip().strip('|').split('|')]
        if len(cells) < 8 or cells[0] in IGNORED_SIGNALS:
            continue
        try:
            max_abs = float(cells[1])
            digits = 16.0 if cells[6] == 'exact' else float(cells[6])
        except ValueError:
            continue
        if worst is None or digits < worst[0]:
            worst = (digits, max_abs)
    return verdict, worst

def build_row(name, status, completion, cases_dir):
    """One table row for one case, or an exit naming the case when its inputs disagree."""
    component, reference = CASE_LABELS[name]
    if status == 'notrun':
        reason = waiver_reason(waiver_path(cases_dir, name))
        return f'| {component} | {reference} | not run | not run | WAIVED, {reason} |'
    verdict = 'PASS' if status == 'passed' else 'FAIL'
    report = os.path.join(cases_dir, name, 'analysis', 'report.md')
    if not os.path.exists(report):
        if verdict == 'PASS':
            sys.exit(f'{name}: passed the run and wrote no report at {report}')
        return f'| {component} | {reference} | {NO_NUMBER} | {NO_NUMBER} | {verdict} |'
    reported, worst = read_report(report)
    if reported != verdict:
        sys.exit(f'{name}: its report says "{reported}" and the run says "{verdict}", '
                 f'so the report is stale')
    if worst is None:
        sys.exit(f'{name}: its report at {report} carries no per-signal statistics')
    return f'| {component} | {reference} | {worst[0]:.1f} | {worst[1]:.2e} | {verdict} |'

def render_region(entries, cases_dir):
    """The region, after policing the label mapping against the run in both directions."""
    registered = {name for name, _, _ in entries}
    for name in sorted(registered - set(CASE_LABELS)):
        print(f'{name}: the run registered it and this script carries no label for it', file=sys.stderr)
    for name in sorted(set(CASE_LABELS) - registered):
        print(f'{name}: this script carries a label for it and the run did not register it', file=sys.stderr)
    if registered != set(CASE_LABELS):
        sys.exit(1)
    return [BEGIN_MARKER, *HEADER, *[build_row(*entry, cases_dir) for entry in entries], END_MARKER]

def replace_region(page, region):
    lines = page.splitlines()
    try:
        begin, end = lines.index(BEGIN_MARKER), lines.index(END_MARKER)
    except ValueError:
        sys.exit('the page carries no generated-region markers, or their text was edited')
    return '\n'.join(lines[:begin] + region + lines[end + 1:]) + '\n'

def check_page(path, regenerated):
    """Regenerate into a temporary copy and compare it against the committed page."""
    with tempfile.TemporaryDirectory() as tmp:
        copy = os.path.join(tmp, os.path.basename(path))
        open(copy, 'w').write(regenerated)
        committed, fresh = open(path).readlines(), open(copy).readlines()
    if committed == fresh:
        return 0
    sys.stderr.writelines(difflib.unified_diff(committed, fresh, f'{path} (committed)',
                                               f'{path} (regenerated)'))
    print(f'{path}: the results region is stale; regenerate it with --write', file=sys.stderr)
    return 1

def self_test():
    """Prove the gate's verdict on every state it can reach, since no case is waived."""
    wrong = 0
    with tempfile.TemporaryDirectory() as root:
        for index, (label, status, completion, waived, expected) in enumerate(GATE_STATES):
            cases_dir = os.path.join(root, str(index), 'cases')
            os.makedirs(os.path.join(cases_dir, 'a_case'))
            if waived:
                open(waiver_path(cases_dir, 'a_case'), 'w').write(
                    'WAIVER_REASON="a synthetic reason"\nWAIVER_DATE="2026-08-06"\n')
            results = os.path.join(cases_dir, 'Test.xml')
            open(results, 'w').write(RESULTS_TEMPLATE.format(name='a_case', status=status,
                                                             completion=completion))
            actual = run_gate(load_entries(results), cases_dir)
            wrong += actual != expected
            print(f'{label:30s} completion="{completion}" expected={expected} '
                  f'actual={actual} {"OK" if actual == expected else "MISMATCH"}')
    return 1 if wrong else 0

def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--gate', action='store_true')
    mode.add_argument('--check', metavar='PAGE')
    mode.add_argument('--write', metavar='PAGE')
    mode.add_argument('--self-test', action='store_true')
    parser.add_argument('--build-dir', default='validation/build')
    parser.add_argument('--validation-dir', default='validation')
    parser.add_argument('--results')
    args = parser.parse_args()

    if args.self_test:
        return self_test()
    cases_dir = os.path.join(args.validation_dir, 'cases')
    entries = load_entries(args.results or find_results(args.build_dir))
    if args.gate:
        return run_gate(entries, cases_dir)

    path = args.check or args.write
    regenerated = replace_region(open(path).read(), render_region(entries, cases_dir))
    if args.check:
        return check_page(path, regenerated)
    open(path, 'w').write(regenerated)
    return 0

if __name__ == '__main__':
    sys.exit(main())
