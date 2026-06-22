#!/usr/bin/env python3
import argparse
import csv
import shutil
import sys

REQUIRED = ["sample", "sample_alias", "condition", "design", "fastq_1", "fastq_2"]


def validate(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        missing = [column for column in REQUIRED if column not in fieldnames]
        if "panel" not in fieldnames and "panel_file" not in fieldnames:
            missing.append("panel or panel_file")
        if missing:
            sys.exit("Missing Pixelator samplesheet columns: " + ", ".join(missing))
        rows = list(reader)
    if not rows:
        sys.exit("Pixelator samplesheet is empty.")
    for row in rows:
        for column in REQUIRED:
            if not row.get(column, "").strip():
                sys.exit(f"Pixelator samplesheet has an empty {column} value.")
        if not row.get("panel", "").strip() and not row.get("panel_file", "").strip():
            sys.exit("Pixelator samplesheet has an empty panel or panel_file value.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    validate(args.input)
    shutil.copyfile(args.input, args.output)


if __name__ == "__main__":
    main()
