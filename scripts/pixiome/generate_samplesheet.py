#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import sys


def csv_value(row, *names, default=""):
    for name in names:
        value = row.get(name, "").strip()
        if value:
            return value
    return default


def all_fastqs(fastq_paths):
    paths = []
    for fastq_path in fastq_paths:
        paths.extend(glob.glob(os.path.join(fastq_path, "**", "*.fastq.gz"), recursive=True))
    return {os.path.abspath(path) for path in paths}


def fastq_matches_sample(path, sample, sample_alias):
    name = os.path.basename(path)
    parts = set(os.path.normpath(path).split(os.sep))
    return (
        name.startswith(f"{sample}_")
        or f"_{sample}_" in name
        or f"Sample_{sample}" in parts
        or (sample_alias and (name.startswith(f"{sample_alias}_") or f"_{sample_alias}_" in name))
    )


def find_fastq_pairs(fastq_set, sample, sample_alias):
    pairs = []
    for r1 in sorted(fastq_set):
        if "_R1_" not in os.path.basename(r1) or not fastq_matches_sample(r1, sample, sample_alias):
            continue
        r2 = r1.replace("_R1_", "_R2_")
        if r2 in fastq_set:
            pairs.append((r1, r2))
    return pairs


def warn_unused_fastqs(fastq_set, used_fastqs):
    unused = sorted(fastq_set - used_fastqs)
    if unused:
        print("WARNING: FASTQ files were found but not included in pixelator_samplesheet.csv:", file=sys.stderr)
        for path in unused:
            print(f"WARNING: {path}", file=sys.stderr)


def write_pixelator_samplesheet_from_libraries(libraries_csv, fastq_paths, output_csv):
    if not os.path.isfile(libraries_csv):
        sys.exit(f"\nlibraries.csv is required for pixiome when --pixelator-samplesheet is not provided: {libraries_csv}\n")
    rows = []
    fastq_set = all_fastqs(fastq_paths)
    used_fastqs = set()
    with open(libraries_csv, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sample = csv_value(row, "sample", "Name")
            sample_alias = csv_value(row, "sample_alias")
            condition = csv_value(row, "condition", "Condition")
            design = csv_value(row, "design")
            panel = csv_value(row, "panel", "Panel")
            panel_file = csv_value(row, "panel_file")
            missing = []
            for column, value in {
                "sample": sample,
                "sample_alias": sample_alias,
                "condition": condition,
                "design": design,
            }.items():
                if not value:
                    missing.append(column)
            if not panel and not panel_file:
                missing.append("panel or panel_file")
            if missing:
                sys.exit("\nlibraries.csv is missing required pixiome values: " + ", ".join(missing) + ".\n")
            if csv_value(row, "fastq_1") and csv_value(row, "fastq_2"):
                pairs = [(os.path.abspath(row["fastq_1"].strip()), os.path.abspath(row["fastq_2"].strip()))]
            else:
                pairs = find_fastq_pairs(fastq_set, sample, sample_alias)
            if not pairs:
                sys.exit(f"\nNo paired FASTQs found for pixiome library '{sample}'.\n")
            for fastq_1, fastq_2 in pairs:
                used_fastqs.update([fastq_1, fastq_2])
                rows.append({
                    "sample": sample,
                    "sample_alias": sample_alias,
                    "condition": condition,
                    "design": design,
                    "panel": panel,
                    "panel_file": os.path.abspath(panel_file) if panel_file else "",
                    "fastq_1": fastq_1,
                    "fastq_2": fastq_2,
                })
    warn_unused_fastqs(fastq_set, used_fastqs)
    with open(output_csv, "w", newline="") as handle:
        panel_columns = []
        if any(row["panel"] for row in rows):
            panel_columns.append("panel")
        if any(row["panel_file"] for row in rows):
            panel_columns.append("panel_file")
        fieldnames = ["sample", "sample_alias", "condition", "design"] + panel_columns + ["fastq_1", "fastq_2"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return sorted({row["sample"] for row in rows})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--libraries", required=True)
    parser.add_argument("--fastq-paths", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    fastq_paths = [path for path in args.fastq_paths.split(",") if path]
    write_pixelator_samplesheet_from_libraries(args.libraries, fastq_paths, args.output)


if __name__ == "__main__":
    main()
