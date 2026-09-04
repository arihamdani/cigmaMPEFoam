#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


NONUNIFORM_PATTERN = re.compile(
    r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
    re.DOTALL,
)
UNIFORM_PATTERN = re.compile(
    r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;"
)


def read_internal_field(path, size=None):
    text = path.read_text()
    match = NONUNIFORM_PATTERN.search(text)

    if match:
        values = [float(value) for value in match.group(2).split()]
        expected_size = int(match.group(1))
        if len(values) != expected_size:
            raise RuntimeError(f"Unexpected field size in {path}")
        return values

    match = UNIFORM_PATTERN.search(text)
    if match and size is not None:
        return [float(match.group(1))] * size

    raise RuntimeError(f"Cannot parse scalar internal field: {path}")


parser = argparse.ArgumentParser()
parser.add_argument("case", type=Path)
parser.add_argument("time")
parser.add_argument("species", nargs="+")
parser.add_argument("--region", default="fluid")
parser.add_argument("--tolerance", type=float, default=1e-7)
args = parser.parse_args()

processor_dirs = sorted(args.case.glob("processor[0-9]*"))
field_roots = processor_dirs if processor_dirs else [args.case]

minimum = float("inf")
maximum = float("-inf")
maximum_location = None

for root in field_roots:
    fields = []
    field_size = None

    for specie in args.species:
        path = root / args.time / args.region / specie
        values = read_internal_field(path, field_size)
        field_size = len(values)
        fields.append(values)

    for celli, values in enumerate(zip(*fields)):
        total = sum(values)
        minimum = min(minimum, total)
        if total > maximum:
            maximum = total
            maximum_location = (root.name, celli, values)

print(f"Internal mass-fraction sum: min={minimum:.12g}, max={maximum:.12g}")
print(
    "Maximum location: "
    f"{maximum_location[0]}, cell={maximum_location[1]}, "
    f"species={maximum_location[2]}"
)

if minimum < 1 - args.tolerance or maximum > 1 + args.tolerance:
    raise SystemExit(
        "Mass-fraction sum is outside "
        f"[1-{args.tolerance:g}, 1+{args.tolerance:g}]"
    )
