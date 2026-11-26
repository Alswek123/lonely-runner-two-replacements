"""Recompute a census range and compare its output with the archived tables."""
import argparse
import csv
from pathlib import Path
import subprocess


def read_expected(filename, first, last):
    with (Path(__file__).resolve().parent / filename).open(newline="") as stream:
        rows = csv.reader(stream)
        next(rows)
        return [row for row in rows if first <= int(row[0]) <= last]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("first", type=int)
    parser.add_argument("last", type=int)
    args = parser.parse_args()
    if not 8 <= args.first <= args.last <= 300:
        parser.error("the archived range is 8 through 300, inclusive")
    command = [str(args.executable.resolve()), str(args.first), str(args.last)]
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    rows = list(csv.reader(result.stdout.splitlines()))
    if any(not row or row[0] not in {"count", "hit"} for row in rows):
        raise ValueError("unexpected census output")
    for tag, filename in [("count", "census_by_n.csv"), ("hit", "tight_sets.csv")]:
        actual = [row[1:] for row in rows if row[0] == tag]
        expected = read_expected(filename, args.first, args.last)
        if actual != expected:
            raise ValueError(f"census differs from {filename}")
    print(f"Archived counts and tight sets agree for n={args.first}..{args.last}.")


if __name__ == "__main__":
    main()
