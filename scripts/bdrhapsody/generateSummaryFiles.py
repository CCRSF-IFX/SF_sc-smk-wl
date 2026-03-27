#!/usr/bin/env python3

import csv
import glob
import os
import re
from pathlib import Path
from shutil import copyfile

import xlsxwriter


METRICS_PATH = "finalreport/"
SUMMARY_PATH = "finalreport/summaries/"
PIPELINE_VERSION_RE = re.compile(r"^(?P<name>.+?) Version (?P<version>[^ ]+)$")
INT_RE = re.compile(r"^[+-]?\d+$")
FLOAT_RE = re.compile(r"^[+-]?(?:\d+\.\d*|\d*\.\d+)$")


def coerce_value(value):
    value = str(value).strip()
    if value in ("", "-"):
        return value, None
    if value.endswith("%"):
        numeric = value[:-1].strip().replace(",", "")
        if INT_RE.match(numeric) or FLOAT_RE.match(numeric):
            return float(numeric) / 100.0, "percent"
        return value, None

    numeric = value.replace(",", "")
    if INT_RE.match(numeric):
        return int(numeric), "int"
    if FLOAT_RE.match(numeric):
        return float(numeric), "float"
    return value, None


def next_nonempty(lines, start_idx):
    idx = start_idx
    while idx < len(lines) and not lines[idx].strip():
        idx += 1
    return idx


def parse_metadata_line(line):
    content = line[2:].strip()
    if not content or set(content) == {"#"}:
        return []
    if " - " in content:
        section, remainder = content.split(" - ", 1)
        entries = []
        for item in remainder.split(" | "):
            if ": " in item:
                key, value = item.split(": ", 1)
                entries.append((f"{section}.{key}", value))
            else:
                entries.append((section, item))
        return entries

    match = PIPELINE_VERSION_RE.match(content)
    if match:
        return [
            ("Pipeline.Name", match.group("name")),
            ("Pipeline.Version", match.group("version")),
        ]

    return [("Pipeline.Info", content)]


def flatten_section(section_name, header, rows):
    flattened = {}
    if not rows:
        return flattened

    id_column = None
    for candidate in ("Library", "Bioproduct_Type"):
        if candidate in header:
            id_column = candidate
            break

    multiple_rows = len(rows) > 1
    for row in rows:
        row_label = None
        if multiple_rows and id_column:
            row_label = row.get(id_column, "")
        for column in header:
            if column == id_column:
                continue
            key = f"{section_name}.{column}"
            if row_label:
                key = f"{section_name}[{row_label}].{column}"
            flattened[key] = row.get(column, "")
    return flattened


def parse_metrics_summary(filename):
    lines = Path(filename).read_text(encoding="utf-8", errors="replace").splitlines()
    metadata = {}
    metrics = {}
    idx = 0

    while idx < len(lines):
        line = lines[idx].strip()
        if not line:
            idx += 1
            continue
        if line.startswith("##"):
            for key, value in parse_metadata_line(line):
                metadata[key] = value
            idx += 1
            continue
        if line.startswith("#"):
            section_name = line.strip("#").strip()
            idx = next_nonempty(lines, idx + 1)
            if idx >= len(lines):
                break
            header_line = lines[idx].strip()
            if not header_line or header_line.startswith("#"):
                continue
            header = next(csv.reader([header_line]))
            idx += 1
            rows = []
            while idx < len(lines):
                candidate = lines[idx].strip()
                if not candidate:
                    idx += 1
                    break
                if candidate.startswith("#"):
                    break
                values = next(csv.reader([candidate]))
                rows.append(dict(zip(header, values)))
                idx += 1
            metrics.update(flatten_section(section_name, header, rows))
            continue
        idx += 1

    return metadata, metrics


def write_sheet(workbook, name, headers, rows, coerce_numbers=True):
    worksheet = workbook.add_worksheet(name)
    format_num = workbook.add_format({"num_format": "#,##0"})
    format_float = workbook.add_format({"num_format": "#,##0.00"})
    format_per = workbook.add_format({"num_format": "0.00%"})
    format_head = workbook.add_format({"bold": True, "italic": True, "text_wrap": True, "align": "center"})

    for col, header in enumerate(headers):
        worksheet.write(0, col, header, format_head)

    for row_idx, row in enumerate(rows, start=1):
        for col_idx, header in enumerate(headers):
            raw_value = row.get(header, "")
            if not coerce_numbers:
                worksheet.write(row_idx, col_idx, raw_value)
                continue
            value, kind = coerce_value(raw_value)
            if kind == "percent":
                worksheet.write(row_idx, col_idx, value, format_per)
            elif kind == "float":
                worksheet.write(row_idx, col_idx, value, format_float)
            elif kind == "int":
                worksheet.write(row_idx, col_idx, value, format_num)
            else:
                worksheet.write(row_idx, col_idx, value)

    worksheet.set_column(0, len(headers) - 1, 22)


def main(output_name="metric_summary"):
    os.makedirs(METRICS_PATH, exist_ok=True)
    os.makedirs(SUMMARY_PATH, exist_ok=True)

    files = sorted(glob.glob("./*/*_Metrics_Summary.csv"))
    workbook = xlsxwriter.Workbook(os.path.join(METRICS_PATH, f"{output_name}.xlsx"))

    summary_headers = ["Sample"]
    metadata_headers = ["Sample"]
    summary_rows = []
    metadata_rows = []

    for filename in files:
        sample = Path(filename).parent.name
        metadata, metrics = parse_metrics_summary(filename)

        summary_row = {"Sample": sample, **metrics}
        metadata_row = {"Sample": sample, **metadata}
        summary_rows.append(summary_row)
        metadata_rows.append(metadata_row)

        for key in metrics:
            if key not in summary_headers:
                summary_headers.append(key)
        for key in metadata:
            if key not in metadata_headers:
                metadata_headers.append(key)

    write_sheet(workbook, "metrics_summary", summary_headers, summary_rows)
    if len(metadata_headers) > 1:
        write_sheet(workbook, "run_metadata", metadata_headers, metadata_rows, coerce_numbers=False)

    workbook.close()

    for filename in sorted(glob.glob("./*/*_Pipeline_Report.html")):
        sample = Path(filename).parent.name
        copyfile(filename, os.path.join(SUMMARY_PATH, f"{sample}_web_summary.html"))


if __name__ == "__main__":
    main()
