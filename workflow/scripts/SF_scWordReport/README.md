# SF_scWordReport (Internal use only)

This repository contains the script to generate a word report for data delivery. 

The script includes the following parts:

1. Use `data/10X_sc.docx` as a template to write a word report in which it includes the project information from the metadata and a run comment to summarize the results from `metric_summary.xlsx`

2. Write an email in word document 

3. Update the tracking information in the csv file below:

`/mnt/ccrsf-ifx/Report_archive/processed_ccrsfifx_SingleCell.csv `
 
## Usage 


```
run_wordreport_sc.py -e data/CS036326/metric_summary.xlsx -m data/CS036326/HN7WMDRX3_Metadata.txt -c /mnt/ccrsf-ifx/Software/tools/GemCode/cellranger-7.1.0/cellranger -p multi -r 240301_A01615_0347_AHN7WMDRX3 --yields 212704.0
```


## TODO list

1. Add unittest coverage for all the code

2. Redesign the format to output the yield information in the tracking sheet. 
