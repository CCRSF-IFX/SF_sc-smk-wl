#!/usr/bin/env python3
"""Build the consolidated final report inputs for the Trekker workflow."""

import argparse
import csv
import shutil
from collections import OrderedDict
from pathlib import Path

import xlsxwriter


def sample_names(libraries):
    """Return unique final sample names in libraries.csv order."""

    with Path(libraries).open(newline="", encoding="utf-8-sig") as handle:
        names = OrderedDict()
        for row in csv.DictReader(handle):
            name = (row.get("Name") or "").strip()
            if name:
                names[name] = None
    if not names:
        raise ValueError(f"No sample names found in {libraries}")
    return list(names)


def read_horizontal_metrics(path):
    """Read Cell Ranger's one-header/one-value metrics CSV."""

    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        try:
            return OrderedDict(next(reader))
        except StopIteration as error:
            raise ValueError(f"Cell Ranger metrics file has no data row: {path}") from error


def read_vertical_metrics(path):
    """Read Trekker's Metrics,Value CSV into one ordered mapping."""

    metrics = OrderedDict()
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["Metrics", "Value"]:
            raise ValueError(
                f"Unexpected Trekker metrics header in {path}: {reader.fieldnames}"
            )
        for row in reader:
            name = (row.get("Metrics") or "").strip()
            if name and name not in {"Sample", "Sample_ID"}:
                metrics[name] = (row.get("Value") or "").strip()
    return metrics


def excel_value(value):
    """Convert formatted metric strings to numbers when conversion is lossless."""

    text = str(value).strip()
    if not text:
        return ""
    numeric = text.replace(",", "")
    if numeric.endswith("%"):
        numeric = numeric[:-1]
    try:
        return float(numeric) if "." in numeric else int(numeric)
    except ValueError:
        return text


def generate_summary(libraries, cellranger_root, trekker_root, output_dir):
    """Merge metrics and copy HTML reports into the standard finalreport layout."""

    cellranger_root = Path(cellranger_root)
    trekker_root = Path(trekker_root)
    output_dir = Path(output_dir)
    summaries_dir = output_dir / "summaries"
    summaries_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    headers = OrderedDict([("Sample", None)])
    for sample in sample_names(libraries):
        cellranger_out = cellranger_root / sample / "outs"
        trekker_out = trekker_root / f"_{sample}" / f"trekker_{sample}" / "output"
        row = OrderedDict([("Sample", sample)])
        row.update(read_horizontal_metrics(cellranger_out / "metrics_summary.csv"))
        row.update(read_vertical_metrics(trekker_out / f"{sample}_summary_metrics.csv"))
        for header in row:
            headers.setdefault(header, None)
        rows.append(row)

        shutil.copy2(
            cellranger_out / "web_summary.html",
            summaries_dir / f"{sample}_web_summary.html",
        )
        shutil.copy2(
            trekker_out / f"{sample}_Trekker_Report.html",
            summaries_dir / f"{sample}_Trekker_Report.html",
        )

    workbook_path = output_dir / "metric_summary.xlsx"
    workbook = xlsxwriter.Workbook(str(workbook_path))
    worksheet = workbook.add_worksheet("metrics_summary")
    header_format = workbook.add_format(
        {"bold": True, "italic": True, "text_wrap": True, "align": "center"}
    )
    integer_format = workbook.add_format({"num_format": "#,##0"})
    decimal_format = workbook.add_format({"num_format": "#,##0.00"})

    header_names = list(headers)
    for column, header in enumerate(header_names):
        worksheet.write(0, column, header, header_format)
        worksheet.set_column(column, column, max(12, min(30, len(header) + 2)))
    for row_index, row in enumerate(rows, start=1):
        for column, header in enumerate(header_names):
            value = excel_value(row.get(header, ""))
            if isinstance(value, int):
                worksheet.write_number(row_index, column, value, integer_format)
            elif isinstance(value, float):
                worksheet.write_number(row_index, column, value, decimal_format)
            else:
                worksheet.write(row_index, column, value)
    workbook.close()
    return workbook_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--libraries", required=True)
    parser.add_argument("--cellranger-root", required=True)
    parser.add_argument("--trekker-root", required=True)
    parser.add_argument("--output-dir", default="finalreport")
    args = parser.parse_args()
    generate_summary(
        args.libraries,
        args.cellranger_root,
        args.trekker_root,
        args.output_dir,
    )


if __name__ == "__main__":
    main()
