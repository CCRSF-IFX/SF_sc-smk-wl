#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

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
    files = glob.glob('./*/all-sample/report/analysis_summary.csv')
    #Filter out aggregate runs if they exist
    #files = [i for i in files if i.split('/')[1] in [j.split('/')[1] for j in glob.glob('./*/*COUNTER*')]]
    files.sort()

    workbook = xlsxwriter.Workbook(metricsPath + arg1+'.xlsx')
    worksheet = workbook.add_worksheet("metrics_summary")
    worksheet.set_column(0, 12, 10.1)
    worksheet.set_column(13, 16, 12.2)
    worksheet.set_column(17,20, 10)

    formatFloat = workbook.add_format({'num_format': '#,##0.00'})
    formatNum = workbook.add_format({'num_format': '#,##'})
    formatPer = workbook.add_format({'num_format': '0.00%'})
    formatHead = workbook.add_format({'bold': True, 'italic': True, 'text_wrap': True, 'align': 'center'})

    # Write headers first
    row = 0
    col = 0
    worksheet.write(row, col, "Sublibrary", formatHead)
    col = 1
    
    # Get statistics from the first file to write column headers
    if files:
        transposed_data = transpose_csv_with_zip(files[0])
        statistics = transposed_data[0]  # First column contains statistic names
        
        for statistic in statistics[1:]:  # Skip first item (which was "statistic")
            worksheet.write(row, col, statistic, formatHead)
            col += 1
    
    row = 1
    sublibraries = list()
    for filename in files:
        print("Processing file: ", filename)
        # Transpose the data using the existing function
        transposed_data = transpose_csv_with_zip(filename)
        
        # First column of transposed data contains statistic names (now as header)
        statistics = transposed_data[0]
        
        # Process each sample (now each row in transposed data)
        for sample_data in transposed_data[1:]:  # Skip first column (statistics)
            sample_name = sample_data[0]  # First value is the sample name
            
            # Process sample name according to your rules
            if sample_name == "combined":
                processed_sample_name = filename.split('/')[1] + "__combined"
            elif sample_name.startswith("all-sample__"):
                processed_sample_name = sample_name.replace("all-sample", filename.split('/')[1])
            else:
                processed_sample_name = sample_name
            print(processed_sample_name)
            # Write sample name in first column
            worksheet.write(row, 0, processed_sample_name)
            
            # Write each statistic value
            col = 1
            for value in sample_data[1:]:  # Skip first value (sample name)
                value = str(value).strip('"')
                try:
                    if '.' in value:
                        worksheet.write(row, col, float(value.replace(',','')), formatFloat)
                    elif '%' in value:
                        worksheet.write(row, col, float(value.strip('%'))/100, formatPer)
                    else:
                        worksheet.write(row, col, int(value.replace(',','')), formatNum)
                except ValueError:
                    # If conversion fails, write as string
                    worksheet.write(row, col, value)
                col += 1
            row += 1

    #for i in sublibraries:
    #    worksheet = workbook.add_worksheet(i)

    workbook.close()

def copyWebSummary():
    try:
        os.makedirs(summaryPath)
    except OSError:
        if not os.path.isdir(summaryPath):
            raise
    files = glob.glob('./*/*_analysis_summary.html')
    for filename in files:
    	copyfile(filename, '%s/%s' % (summaryPath, filename.split('/')[1] + "_" + filename.split('/')[-1]))

# Example functions for transposing CSV data after csv.reader

def transpose_csv_with_zip(filename):
    """Transpose entire CSV using zip(*data)"""
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        data = list(reader)  # Read all rows into a list
        transposed = list(zip(*data))  # Transpose: rows become columns
        return transposed

if __name__ == "__main__":
    if len(sys.argv) == 1:
        main()
    else:
        main(sys.argv[1])
