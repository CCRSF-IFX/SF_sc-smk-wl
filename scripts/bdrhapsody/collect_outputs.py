#!/usr/bin/env python3

import argparse
import csv
import gzip
import html
import shutil
import sys
import zipfile
from pathlib import Path


def find_one(rawdir, pattern):
    matches = sorted(rawdir.rglob(pattern))
    return matches[0] if matches else None


def sum_stats(bioproduct_stats):
    metrics = {}
    with bioproduct_stats.open(newline="", encoding="utf-8") as handle:
        rows = [line for line in handle if not line.startswith("#")]
    reader = csv.DictReader(rows)
    totals = {
        "Raw_Reads": 0,
        "Raw_Molecules": 0,
        "RSEC_Adjusted_Molecules": 0,
        "RSEC_Adjusted_Reads_non-singleton": 0,
        "RSEC_Adjusted_Molecules_non-singleton": 0,
        "DBEC_Adjusted_Reads": 0,
        "DBEC_Adjusted_Molecules": 0,
    }
    for row in reader:
        for key in totals:
            totals[key] += int(float(row.get(key, 0) or 0))
    metrics.update(totals)
    return metrics


def count_gzip_lines(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return sum(1 for _ in handle)


def render_web_summary(sample, metrics, rawdir):
    rows = "\n".join(
        f"<tr><th>{html.escape(key)}</th><td>{html.escape(str(value))}</td></tr>"
        for key, value in metrics.items()
    )
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(sample)} BD Rhapsody Summary</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; }}
    table {{ border-collapse: collapse; }}
    th, td {{ border: 1px solid #ccc; padding: 0.5rem 0.75rem; text-align: left; }}
    th {{ background: #f5f5f5; }}
    code {{ background: #f3f3f3; padding: 0.1rem 0.3rem; }}
  </style>
</head>
<body>
  <h1>{html.escape(sample)} BD Rhapsody Summary</h1>
  <p>Raw vendor output directory: <code>{html.escape(str(rawdir))}</code></p>
  <table>
    <tbody>
      {rows}
    </tbody>
  </table>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--rawdir", required=True)
    parser.add_argument("--outs", required=True)
    parser.add_argument("--vendor-log", required=True)
    parser.add_argument("--vendor-err", required=True)
    args = parser.parse_args()

    rawdir = Path(args.rawdir)
    outs = Path(args.outs)
    matrix_dir = outs / "filtered_feature_bc_matrix"
    vendor_dir = outs / "vendor"
    matrix_dir.mkdir(parents=True, exist_ok=True)
    vendor_dir.mkdir(parents=True, exist_ok=True)

    rsec_zip = find_one(rawdir, "*_RSEC_MolsPerCell_MEX.zip")
    bioproduct_stats = find_one(rawdir, "*_Bioproduct_Stats.csv")
    if rsec_zip is None or bioproduct_stats is None:
        sys.stderr.write(
            "BD Rhapsody outputs are incomplete. "
            "Expected both *_RSEC_MolsPerCell_MEX.zip and *_Bioproduct_Stats.csv.\n"
        )
        sys.exit(1)

    with zipfile.ZipFile(rsec_zip) as zf:
        zf.extractall(matrix_dir)

    for path in rawdir.iterdir():
        if path.is_file():
            shutil.copy2(path, vendor_dir / path.name)

    shutil.copy2(args.vendor_log, vendor_dir / Path(args.vendor_log).name)
    shutil.copy2(args.vendor_err, vendor_dir / Path(args.vendor_err).name)

    metrics = sum_stats(bioproduct_stats)
    metrics["Putative_Cells"] = count_gzip_lines(matrix_dir / "barcodes.tsv.gz")
    metrics["Features"] = count_gzip_lines(matrix_dir / "features.tsv.gz")

    metrics_path = outs / "metrics_summary.csv"
    with metrics_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics.keys()))
        writer.writeheader()
        writer.writerow(metrics)

    web_summary_path = outs / "web_summary.html"
    web_summary_path.write_text(
        render_web_summary(args.sample, metrics, rawdir),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
