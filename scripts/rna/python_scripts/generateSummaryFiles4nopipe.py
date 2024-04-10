#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Feb 23 2022

@author: Shaojun Xie

metric_summary is used because previous developer used this name for 10x data
usage: python generateSummaryFiles.py metric_summary

"""


import glob, xlsxwriter, csv, sys, ntpath, os
import re
from shutil import copyfile

metricsPath = 'finalreport/'
summaryPath = 'finalreport/summaries/'
def main(arg1='metric_summary'):
    createMetricsSummary(arg1)

def createMetricsSummary(arg1):
    try:
        os.makedirs(metricsPath)
    except OSError:
        if not os.path.isdir(metricsPath):
            raise

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
    with open('multiqc4report_data/multiqc_general_stats.txt', 'r') as csvfile:
        f = csv.reader(csvfile, delimiter='\t', quotechar='"')
        header = next(f)
        worksheet.write(0, 0, header[0])          
        worksheet.write(0, 1, re.sub("FastQC_mqc-generalstats-fastqc-total_sequences", "FastQC_total_sequences", header[-1])) 
        for line in f: 
            worksheet.write(row, 0, line[0]) 
            worksheet.write(row, 1, int(float(line[-1])), formatNum)
            row += 1
    workbook.close()

if __name__ == "__main__":
    if len(sys.argv) == 1:
        main()
    else:
        main(sys.argv[1])
