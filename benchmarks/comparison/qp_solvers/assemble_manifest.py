#!/usr/bin/env python3
"""Assemble a run directory's manifest.json from the pre-registered template and the harness's outputs.

The template (manifest_template.json, a verbatim copy of the pre-registration's manifest template) is
filled where it holds null and is otherwise left as it is, except that each arm's tier object gains the
arm's remaining effective settings as extra fields. Before anything is written, the harness's own
description of what it runs (settings, constants, expressibility) is checked against the template; a
disagreement stops the assembly, so a manifest is only ever written for a harness that runs what was
pre-registered.

Standard-library Python 3 only.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import struct
import sys
from pathlib import Path

CSV_NAME = re.compile(r"^(dense|sparse)_rep([0-9]+)\.csv\.gz$")
PROBE_NAME = re.compile(r"^(dense|sparse)_rep([0-9]+)_probes\.json$")


class AssemblyError(Exception):
    pass


def same(a: object, b: object) -> bool:
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(map(same, a, b))
    return type(a) is type(b) and a == b


def load(path: Path) -> object:
    return json.loads(path.read_text(encoding="utf-8"))


def merged_settings(template: dict, described: dict) -> dict:
    if template.keys() != described.keys():
        raise AssemblyError(f"arm_settings: the harness runs arms {sorted(described)}, the template lists {sorted(template)}")
    merged = {}
    for arm, tiers in template.items():
        if tiers.keys() != described[arm].keys():
            raise AssemblyError(f"arm_settings: {arm} runs tiers {sorted(described[arm])}, the template lists {sorted(tiers)}")
        merged[arm] = {}
        for tier, fields in tiers.items():
            effective = described[arm][tier]
            for field, value in fields.items():
                if not same(effective.get(field), value):
                    raise AssemblyError(f"arm_settings: {arm} {tier} {field} is {effective.get(field)!r}, pre-registered {value!r}")
            merged[arm][tier] = {**effective, **fields}
    return merged


def check_constants(template: dict, described: dict) -> None:
    if not same(template, described):
        raise AssemblyError("constants: the harness's constants differ from the pre-registered ones\n"
                            f"  harness:  {json.dumps(described, sort_keys=True)}\n"
                            f"  template: {json.dumps(template, sort_keys=True)}")


def check_expressibility(template: list, described: dict) -> None:
    harness = {(cohort_arm.split("/")[0], cohort_arm.split("/")[1], cell, flag)
               for cohort_arm, cells in described.items() for cell, flag in cells.items()}
    registered = {(e["cohort"], e["arm"], e["cell"], e["expressible"]) for e in template}
    if harness != registered:
        raise AssemblyError(f"expressibility: differs from the template in {sorted(harness ^ registered)}")


def terminal_weight_hashes(described: dict) -> dict:
    hashes = {}
    for cell, entries in described.items():
        values = [entries[str(i)] for i in range(len(entries) - 1)]
        hashes[cell] = hashlib.sha256(struct.pack(f"<{len(values)}d", *values)).hexdigest()
    return hashes


def merged_probes(run_dir: Path) -> dict:
    probes: dict = {}
    files = sorted(p for p in run_dir.iterdir() if PROBE_NAME.match(p.name))
    if not files:
        raise AssemblyError(f"{run_dir}: no probe records")
    for path in files:
        for arm, result in load(path)["probes"].items():
            error = result["max_dual_error"]
            error = float("inf") if error is None else float(error)
            previous = probes.get(arm, {"passed": True, "max_dual_error": 0.0})
            probes[arm] = {"passed": previous["passed"] and result["passed"],
                           "max_dual_error": max(previous["max_dual_error"], error)}
    for arm, result in probes.items():
        if result["max_dual_error"] == float("inf"):
            result["max_dual_error"] = sys.float_info.max
    return probes


def csv_hashes(run_dir: Path) -> dict:
    files = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(run_dir.iterdir()) if CSV_NAME.match(p.name)}
    if not files:
        raise AssemblyError(f"{run_dir}: no <cohort>_rep<r>.csv.gz files")
    for name in files:
        with gzip.open(run_dir / name, "rt", encoding="utf-8") as stream:
            stream.readline()
    return files


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--describe", type=Path, required=True)
    parser.add_argument("--gate", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--fill", type=Path, required=True,
                        help="JSON object with run_id, started_utc, finished_utc, preregistration_commit, argmin_sha, "
                             "competitor_pins, ctrlpp_harness_sha, ctrlpp_harness_dirty, compiler_version, c_flags, "
                             "cxx_flags, host and clock_state")
    return parser.parse_args()


def assemble(args: argparse.Namespace) -> dict:
    manifest = load(args.template)
    described = load(args.describe)
    fill = load(args.fill)
    check_constants(manifest["constants"], described["constants"])
    check_expressibility(manifest["expressibility"], described["expressible"])
    manifest["arm_settings"] = merged_settings(manifest["arm_settings"], described["arm_settings"])
    manifest["update_api"] = described["update_api"]
    manifest["terminal_weights"] = terminal_weight_hashes(described["terminal_weights"])
    manifest["binding_gate"] = load(args.gate)
    manifest["probes"] = merged_probes(args.run_dir)
    manifest["csv_files"] = csv_hashes(args.run_dir)
    for key in ("run_id", "started_utc", "finished_utc", "preregistration_commit", "argmin_sha", "competitor_pins",
                "ctrlpp_harness_sha", "ctrlpp_harness_dirty", "host", "clock_state"):
        manifest[key] = fill[key]
    for key in ("compiler_version", "c_flags", "cxx_flags"):
        manifest["build"][key] = fill[key]
    return manifest


def main() -> int:
    args = arguments()
    try:
        manifest = assemble(args)
    except (AssemblyError, KeyError, OSError, ValueError) as error:
        print(f"assemble_manifest: {error}", file=sys.stderr)
        return 1
    missing = [key for key, value in manifest.items() if value is None]
    if missing:
        print(f"assemble_manifest: still null: {', '.join(missing)}", file=sys.stderr)
        return 1
    text = json.dumps(manifest, sort_keys=True, indent=2, allow_nan=False) + "\n"
    (args.run_dir / "manifest.json").write_text(text, encoding="utf-8", newline="\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
