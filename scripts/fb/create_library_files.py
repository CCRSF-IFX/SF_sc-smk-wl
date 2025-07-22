#!/usr/bin/env /is2/projects/CCR-SF/active/Software/tools/Anaconda/3.6/install/bin/python

# /mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/single_cell/featurebarcode/scripts/python_scripts/create_library_files.py 

import argparse, sys

class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

from argparse import ArgumentParser

def main(raw_args=None):
    parser = argparse.ArgumentParser(description="Create libraries file for each sample given input csv", formatter_class=SmartFormatter)
    parser.add_argument('file_name', metavar='libraries.csv', action = "store",
        type=str, help="CSV file containing the sample name, fastq path, and the library type for each file. Expected column order in file: \n(Sample) Name, Flowcell, Sample, Library Type")
    parser.add_argument('fastqs', metavar="outs/fastq_path/HLGW3DRXX",
        nargs='?', action="store", type=str,
        help="Full path to FASTQ files, used to help fill in new library files. If missing will just use values in the libraries.csv file. Multiple paths should be comma delimited")
    # parser.add_argument('--ocm', action='store_true', default=False,
    #     help="Use this flag if OCM platform is used")

    args = parser.parse_args(raw_args)
    if args.fastqs != None:
        fastqs = args.fastqs.split(',')

    ## Check if the input file is in OCM format
    import pandas as pd
    df = pd.read_csv(args.file_name, header=0)
    if 'ocm_barcode_ids' in df.columns:
        if df['ocm_barcode_ids'].any():
            print("OCM platform detected, using OCM specific libraries file format")
            args.ocm = True
        else:
            args.ocm = False
    if args.ocm == False: 
        with open(args.file_name) as f:
            headers = next(f).strip().split(',')
            #print(headers)
            samples = dict()
            for line in f:
                line = line.strip().split(',')
                if line[0] in samples:
                    samples[line[0]].append(line[1:])
                else:
                    samples[line[0]] = [line[1:]]

        for sample in samples:
            text = []
            for values in samples[sample]:
                if args.fastqs != None:
                    runs = [path for path in fastqs if values[0] in path]
                    if len(runs) != 1:
                        sys.exit("Problems finding unique match for %s in %s" % (values[0], args.fastqs))
                    else:
                        if values[2] != values[-1]:
                            text.append(",".join([runs[0], values[1], values[2], values[-1]]))
                        else:
                            text.append(",".join([runs[0], values[1], values[2]]))
                else:
                    text.append(",".join(values))

            with open('%s_libraries.csv' % sample, 'w') as f:
                f.write('fastqs,sample,library_type\n')
                f.write('\n'.join(text))
    else:
        import pandas as pd
        df = pd.read_csv(args.file_name, header=0)
        samples = df['Sample'].unique()
        for sample in samples:
            subdf = df[df['Sample'] == sample].copy()
            # If fastqs is provided, try to match each row's Flowcell to a fastq path
            if args.fastqs is not None:
                matched_fastqs = []
                for idx, row in subdf.iterrows():
                    # Try to match Flowcell in fastqs path
                    matches = [fq for fq in fastqs if str(row['Flowcell']) in fq]
                    if len(matches) == 1:
                        matched_fastqs.append(matches[0])
                    elif len(matches) > 1:
                        sys.exit(f"Multiple matches for Flowcell '{row['Flowcell']}' in fastqs: {matches}")
                    else:
                        matched_fastqs.append("")  # Or handle as needed
                subdf.insert(0, 'fastqs', matched_fastqs)
            if "Name" in subdf.columns:
                subdf = subdf.rename(columns={"Name": "sample"})
            if "Type" in subdf.columns:
                subdf = subdf.rename(columns={"Type": "library_type"})
            # Remove "Flowcell" and "Sample" columns if present
            for col in ["Flowcell", "Sample"]:
                if col in subdf.columns:
                    subdf = subdf.drop(columns=[col])
            # Write all columns for this sample to a new CSV
            out_file = f"{sample}_libraries.csv"
            subdf.to_csv(out_file, index=False)
    #print(samples)


if __name__ == '__main__':
    main()
