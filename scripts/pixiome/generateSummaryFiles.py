#!/usr/bin/env python3
import glob
import html
import json
import os
import re
from pathlib import Path
from shutil import copyfile

import xlsxwriter

METRICS_PATH = "finalreport"
SUMMARY_PATH = os.path.join(METRICS_PATH, "summaries")
SUMMARY_HTML = "pixelator/experiment-summary.html"
SAMPLE_METRICS_SIGNATURE = ("Sample ID", "Number of cells in sample")


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


def clean_html_text(value):
    text = re.sub(r"<i\b.*?</i>", "", value, flags=re.S)
    text = re.sub(r"<[^>]+>", "", text)
    return " ".join(html.unescape(text).split())


def transpose_datatable(data):
    if not data:
        return []
    row_count = max(len(column) for column in data)
    rows = []
    for row_idx in range(row_count):
        row = []
        for column in data:
            row.append(column[row_idx] if row_idx < len(column) else "")
        rows.append(row)
    return rows


def normalize_cell(value):
    if isinstance(value, (dict, list)):
        return json.dumps(value, sort_keys=True)
    return value


def sample_metrics_from_html(path=SUMMARY_HTML):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Pixelator HTML summary not found: {path}")

    text = Path(path).read_text()
    pattern = re.compile(r'<script type="application/json" data-for="[^"]+">(.*?)</script>', re.S)
    for match in pattern.finditer(text):
        try:
            payload = json.loads(html.unescape(match.group(1)))
        except json.JSONDecodeError:
            continue

        table = payload.get("x", {})
        if "container" not in table or "data" not in table:
            continue

        headers = [clean_html_text(header) for header in re.findall(r"<th>(.*?)</th>", table["container"], re.S)]
        if tuple(headers[: len(SAMPLE_METRICS_SIGNATURE)]) != SAMPLE_METRICS_SIGNATURE:
            continue

        rows = [[normalize_cell(value) for value in row] for row in transpose_datatable(table["data"])]
        return headers, rows

    raise ValueError("Could not find the Pixelator Sample metrics table in experiment-summary.html")


def write_sheet(workbook, name, headers, rows):
    worksheet = workbook.add_worksheet(name)
    head = workbook.add_format({"bold": True, "italic": True, "text_wrap": True, "align": "center"})
    for col, header in enumerate(headers):
        worksheet.write(0, col, header, head)
    for row_idx, row in enumerate(rows, start=1):
        if isinstance(row, dict):
            values = [row.get(header, "") for header in headers]
        else:
            values = row
        for col_idx, value in enumerate(values):
            worksheet.write(row_idx, col_idx, value)
    worksheet.set_column(0, len(headers) - 1, 24)


def main(output_name="metric_summary"):
    os.makedirs(SUMMARY_PATH, exist_ok=True)
    headers, rows = sample_metrics_from_html()

    workbook = xlsxwriter.Workbook(os.path.join(METRICS_PATH, f"{output_name}.xlsx"))
    write_sheet(workbook, "metrics_summary", headers, rows)
    write_sheet(
        workbook,
        "run_metadata",
        ["Pipeline.Name", "Pipeline.Version", "SampleSheet"],
        [{"Pipeline.Name": "Pixelator", "Pipeline.Version": pixelator_version(), "SampleSheet": os.path.abspath("samplesheet.pixiome.csv")}],
    )
    workbook.close()

    if os.path.exists(SUMMARY_HTML):
        copyfile(SUMMARY_HTML, os.path.join(SUMMARY_PATH, "experiment-summary.html"))


if __name__ == "__main__":
    main()
