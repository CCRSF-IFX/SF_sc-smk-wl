#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Feb 2, 2023

@author: Shaojun Xie

usage: python generateSummaryFiles.py metric_summary

"""


import glob, xlsxwriter, csv, sys, ntpath, os
from shutil import copyfile

metricsPath = 'finalreport/'
summaryPath = 'finalreport/summaries/'
def main(arg1='metric_summary'):
    createMetricsSummary(arg1)
    copyWebSummary()

def createMetricsSummary(arg1):
    try:
        os.makedirs(metricsPath)
    except OSError:
        if not os.path.isdir(metricsPath):
            raise
    files = glob.glob('./*/metrics/*/metrics_summary.csv')
    #Filter out aggregate runs if they exist
    #files = [i for i in files if i.split('/')[1] in [j.split('/')[1] for j in glob.glob('./*/*COUNTER*')]]
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
    for filename in files:
        with open(filename, 'r') as csvfile:
            f = csv.reader(csvfile, delimiter=',', quotechar='"')
            header = next(f)
            line = next(f)
            worksheet.write(row, 0, filename.split('/')[1])
            samples.append(filename.split('/')[1])
            worksheet.write(row, 1, filename.split('/')[3]) 
            col = 2
            for i in line:
                i = i.strip('"')
                if '%' in i or '.' in i:
                    worksheet.write(row, col, float(i.strip('%'))/100, formatPer)
                else:
                    worksheet.write(row, col, int(i.replace(',','')), formatNum)
                col += 1
            row += 1

    col = 2
    row = 0
    worksheet.write(0, 0, "sample", formatHead)
    worksheet.write(0, 1, "sensitivity_level", formatHead)
    for i in header:
        worksheet.write(row, col, i, formatHead)
        col += 1

    #for i in samples:
    #    worksheet = workbook.add_worksheet(i)

    workbook.close()

def copyWebSummary():
    try:
        os.makedirs(summaryPath)
    except OSError:
        if not os.path.isdir(summaryPath):
            raise
    files = glob.glob('./*/report.html')
    for filename in files:
        basename_html = os.path.basename(filename) 
    	copyfile(filename, '%s/%s_%s' % (summaryPath, filename.split('/')[1], basename_html))

if __name__ == "__main__":
    if len(sys.argv) == 1:
        main()
    else:
        main(sys.argv[1])