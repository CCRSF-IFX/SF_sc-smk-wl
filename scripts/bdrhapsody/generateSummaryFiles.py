#!/usr/bin/env python3

import csv
import glob
import os
from shutil import copyfile

import xlsxwriter


METRICS_PATH = "finalreport/"
SUMMARY_PATH = "finalreport/summaries/"


def coerce_value(value):
    value = str(value).strip()
    if value == "":
        return value, None
    if value.endswith("%"):
        try:
            return float(value[:-1]) / 100.0, "percent"
        except ValueError:
            return value, None
    try:
        if "." in value:
            return float(value), "float"
        return int(value), "int"
    except ValueError:
        return value, None


def main(output_name="metric_summary"):
    os.makedirs(METRICS_PATH, exist_ok=True)
    os.makedirs(SUMMARY_PATH, exist_ok=True)

    files = sorted(glob.glob("./*/outs/metrics_summary.csv"))
    workbook = xlsxwriter.Workbook(os.path.join(METRICS_PATH, f"{output_name}.xlsx"))
    worksheet = workbook.add_worksheet("metrics_summary")

    format_num = workbook.add_format({"num_format": "#,##0"})
    format_float = workbook.add_format({"num_format": "#,##0.00"})
    format_per = workbook.add_format({"num_format": "0.00%"})
    format_head = workbook.add_format({"bold": True, "italic": True, "text_wrap": True, "align": "center"})

    headers = ["Sample"]
    rows = []
    for filename in files:
        with open(filename, newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            data = next(reader)
        sample = filename.split("/")[1]
        for key in data:
            if key not in headers:
                headers.append(key)
        rows.append((sample, data))

    for col, header in enumerate(headers):
        worksheet.write(0, col, header, format_head)

    for row_idx, (sample, data) in enumerate(rows, start=1):
        worksheet.write(row_idx, 0, sample)
        for col_idx, header in enumerate(headers[1:], start=1):
            value, kind = coerce_value(data.get(header, ""))
            if kind == "percent":
                worksheet.write(row_idx, col_idx, value, format_per)
            elif kind == "float":
                worksheet.write(row_idx, col_idx, value, format_float)
            elif kind == "int":
                worksheet.write(row_idx, col_idx, value, format_num)
            else:
                worksheet.write(row_idx, col_idx, value)

    worksheet.set_column(0, len(headers) - 1, 18)
    workbook.close()

    for filename in sorted(glob.glob("./*/outs/web_summary.html")):
        sample = filename.split("/")[1]
        copyfile(filename, os.path.join(SUMMARY_PATH, f"{sample}_web_summary.html"))


if __name__ == "__main__":
    main()
