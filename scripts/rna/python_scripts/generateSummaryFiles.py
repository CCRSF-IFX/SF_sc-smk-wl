#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Jul 31 2018

@author: Vicky

usage: python generateSummaryFiles.py metric_summary

"""


import glob, xlsxwriter, csv, sys, ntpath, os
from shutil import copyfile

metricsPath = 'finalreport/'
summaryPath = 'finalreport/summaries/'
def main(arg1='metric_summary'):
    workbook = createMetricsSummary(arg1)
    copyWebSummary()
    create_cell_type_summary(workbook)
    workbook.close()

def create_cell_type_summary(workbook):
    import collections
    files = glob.glob("./*/outs/cell_types/cell_types.csv")
    cell_type_counts = collections.defaultdict(lambda: collections.defaultdict(int))
    fine_type_counts = collections.defaultdict(lambda: collections.defaultdict(int))
    for filepath in files:
        sample = filepath.split('/')[1]
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                coarse = row['coarse_cell_type']
                fine = row['fine_cell_type']
                count = int(row.get('cell_count_in_model', 1))
                cell_type_counts[sample][coarse] += 1
                fine_type_counts[sample][fine] += 1
    # Calculate total counts for sorting
    coarse_totals = collections.Counter()
    for sample in cell_type_counts:
        for ct, val in cell_type_counts[sample].items():
            coarse_totals[ct] += val
    sorted_coarse_types = [ct for ct, _ in coarse_totals.most_common()]
    fine_totals = collections.Counter()
    for sample in fine_type_counts:
        for ft, val in fine_type_counts[sample].items():
            fine_totals[ft] += val
    sorted_fine_types = [ft for ft, _ in fine_totals.most_common()]
    # Sort sample names alphabetically
    sorted_samples = sorted(cell_type_counts.keys())
    sorted_fine_samples = sorted(fine_type_counts.keys())
    # Add new sheet for coarse cell types
    worksheet = workbook.add_worksheet("cell_type_summary")
    worksheet.write(0, 0, "Sample")
    for col, ct in enumerate(sorted_coarse_types, 1):
        worksheet.write(0, col, ct)
    for row, sample in enumerate(sorted_samples, 1):
        worksheet.write(row, 0, sample)
        for col, ct in enumerate(sorted_coarse_types, 1):
            worksheet.write(row, col, cell_type_counts[sample].get(ct, 0))
    # Add new sheet for fine cell types
    worksheet2 = workbook.add_worksheet("fine_cell_type_summary")
    worksheet2.write(0, 0, "Sample")
    for col, ft in enumerate(sorted_fine_types, 1):
        worksheet2.write(0, col, ft)
    for row, sample in enumerate(sorted_fine_samples, 1):
        worksheet2.write(row, 0, sample)
        for col, ft in enumerate(sorted_fine_types, 1):
            worksheet2.write(row, col, fine_type_counts[sample].get(ft, 0))

def createMetricsSummary(arg1):
    try:
        os.makedirs(metricsPath)
    except OSError:
        if not os.path.isdir(metricsPath):
            raise
    files = glob.glob('./*/outs/metrics_summary.csv')
    #Filter out aggregate runs if they exist
    files = [i for i in files if i.split('/')[1] in [j.split('/')[1] for j in glob.glob('./*/*COUNTER*')]]
    files.sort()

    workbook = xlsxwriter.Workbook(metricsPath + arg1+'.xlsx')
    worksheet = workbook.add_worksheet("metrics_summary")
    worksheet.set_column(0, 12, 10.1)
    worksheet.set_column(13, 16, 12.2)
    worksheet.set_column(17,20, 10)

    formatNum = workbook.add_format({'num_format': '#,###'})
    formatPer = workbook.add_format({'num_format': '0.00%'})
    formatHead = workbook.add_format({'bold': True, 'italic': True, 'text_wrap': True, 'align': 'center'})

    row = 1
    samples = list()
    header = ""
    for filename in files:
        with open(filename, 'r') as csvfile:
            f = csv.reader(csvfile, delimiter=',', quotechar='"')
            header = next(f)
            line = next(f)
            worksheet.write(row, 0, filename.split('/')[1])
            samples.append(filename.split('/')[1])
            col = 1
            for i in line:
                i = i.strip('"')
                if '%' in i:
                    worksheet.write(row, col, float(i.strip('%'))/100, formatPer)
                else:
                    worksheet.write(row, col, int(i.replace(',','')), formatNum)
                col += 1
            row += 1

    col = 1
    row = 0
    worksheet.write(0, 0, "Sample", formatHead)
    for i in header:
        worksheet.write(row, col, i, formatHead)
        col += 1

    #for i in samples:
    #    worksheet = workbook.add_worksheet(i)

    # Do NOT close the workbook here, return the open workbook object
    return workbook

def copyWebSummary():
    try:
        os.makedirs(summaryPath)
    except OSError:
        if not os.path.isdir(summaryPath):
            raise
    files = glob.glob('./*/outs/web_summary.html')
    for filename in files:
    	copyfile(filename, '%s/%s_web_summary.html' % (summaryPath, filename.split('/')[1]))
    files_cell_type = glob.glob('./*/outs/cell_types/web_summary_cell_types.html')
    for filename in files_cell_type:
    	copyfile(filename, '%s/%s_web_summary_cell_types.html' % (summaryPath, filename.split('/')[1]))

if __name__ == "__main__":
    if len(sys.argv) == 1:
        main()
    else:
        main(sys.argv[1])
