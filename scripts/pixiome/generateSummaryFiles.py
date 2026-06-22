#!/usr/bin/env python3
import csv
import glob
import json
import os
import re
from pathlib import Path
from shutil import copyfile

import xlsxwriter

METRICS_PATH = "finalreport"
SUMMARY_PATH = os.path.join(METRICS_PATH, "summaries")


def flatten(prefix, value, out):
    if isinstance(value, dict):
        for key, nested in value.items():
            flatten(f"{prefix}.{key}" if prefix else key, nested, out)
    else:
        out[prefix] = value


def samplesheet_rows(path="samplesheet.pixiome.csv"):
    rows = {}
    if not os.path.exists(path):
        return rows
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            rows.setdefault(row["sample"], row)
    return rows


def pixelator_version():
    container = ""
    for params_file in sorted(glob.glob("pipeline_info/params_*.json")):
        with open(params_file) as handle:
            params = json.load(handle)
        container = params.get("pixelator_container", container)
    if not container and os.path.exists("params.pixiome.yaml"):
        text = Path("params.pixiome.yaml").read_text()
        match = re.search(r"pixelator_container:\s*[\"']?([^\"'\n]+)", text)
        if match:
            container = match.group(1)
    match = re.search(r":(\d+\.\d+\.\d+)", container)
    return match.group(1) if match else "UNKNOWN"


def write_sheet(workbook, name, headers, rows):
    worksheet = workbook.add_worksheet(name)
    head = workbook.add_format({"bold": True, "italic": True, "text_wrap": True, "align": "center"})
    num = workbook.add_format({"num_format": "#,##0.0000"})
    for col, header in enumerate(headers):
        worksheet.write(0, col, header, head)
    for row_idx, row in enumerate(rows, start=1):
        for col_idx, header in enumerate(headers):
            value = row.get(header, "")
            if isinstance(value, float):
                worksheet.write(row_idx, col_idx, value, num)
            else:
                worksheet.write(row_idx, col_idx, value)
    worksheet.set_column(0, len(headers) - 1, 22)


def main(output_name="metric_summary"):
    os.makedirs(SUMMARY_PATH, exist_ok=True)
    sheet_rows = samplesheet_rows()
    reports = sorted(glob.glob("pixelator/analysis/*.report.json"))
    headers = [
        "Sample",
        "sample_alias",
        "condition",
        "design",
        "panel",
        "product_id",
        "report_type",
        "k_cores.median_average_k_core",
        "svd.median_variance_explained_3d",
    ]
    rows = []
    for report in reports:
        with open(report) as handle:
            data = json.load(handle)
        flat = {}
        flatten("", data, flat)
        sample = flat.get("sample_id", Path(report).name.split(".")[0])
        row = {"Sample": sample}
        for column in ["sample_alias", "condition", "design", "panel"]:
            row[column] = sheet_rows.get(sample, {}).get(column, "")
        row.update(flat)
        rows.append(row)
        for key in flat:
            if key not in headers and key != "sample_id":
                headers.append(key)

    workbook = xlsxwriter.Workbook(os.path.join(METRICS_PATH, f"{output_name}.xlsx"))
    write_sheet(workbook, "metrics_summary", headers, rows)
    write_sheet(
        workbook,
        "run_metadata",
        ["Pipeline.Name", "Pipeline.Version", "SampleSheet"],
        [{"Pipeline.Name": "Pixelator", "Pipeline.Version": pixelator_version(), "SampleSheet": os.path.abspath("samplesheet.pixiome.csv")}],
    )
    workbook.close()

    if os.path.exists("pixelator/experiment-summary.html"):
        copyfile("pixelator/experiment-summary.html", os.path.join(SUMMARY_PATH, "experiment-summary.html"))


if __name__ == "__main__":
    main()
