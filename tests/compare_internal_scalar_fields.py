#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


FIELD_PATTERN = re.compile(
    r"internalField\s+nonuniform\s+List<scalar>\s+\d+\s*\((.*?)\)\s*;",
    re.DOTALL,
)


def read_internal_field(path):
    match = FIELD_PATTERN.search(Path(path).read_text())
    if not match:
        raise RuntimeError(f"Cannot parse nonuniform scalar field: {path}")
    return [float(value) for value in match.group(1).split()]


parser = argparse.ArgumentParser()
parser.add_argument("reference")
parser.add_argument("candidate")
parser.add_argument("fields", nargs="+")
parser.add_argument("--relative-tolerance", type=float, default=1e-4)
parser.add_argument("--absolute-tolerance", type=float, default=1e-12)
args = parser.parse_args()

for field in args.fields:
    reference = read_internal_field(Path(args.reference) / field)
    candidate = read_internal_field(Path(args.candidate) / field)

    if len(reference) != len(candidate):
        raise RuntimeError(f"Field size differs for {field}")

    max_difference = max(abs(a - b) for a, b in zip(reference, candidate))
    scale = max(abs(value) for value in reference)
    tolerance = args.absolute_tolerance + args.relative_tolerance * scale

    if max_difference > tolerance:
        raise RuntimeError(
            f"{field}: max difference {max_difference} exceeds {tolerance}"
        )

    print(f"{field}: max difference = {max_difference:.6e}")

