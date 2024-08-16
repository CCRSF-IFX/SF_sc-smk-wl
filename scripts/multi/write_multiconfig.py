#!/usr/bin/env python

import argparse, csv, re

def main(raw_args=None):
    parser = argparse.ArgumentParser(description="""Help to set up and run the single cell multi pipeline""")
    parser.add_argument("-o", "--output", metavar="output.csv",
        action = "store", type=str, required=True,
        help="Output file name")
    parser.add_argument("-r", "--ref", metavar="refdata-gex-GRCh38-2020-A",
        action = "store", type=str, required=True,
        help="Path to reference")
    parser.add_argument("-l", "--lib", metavar="libraries.csv",
        action = "store", type=str, required=True,
        help="Path to libraries file to create the config file based off of")
    parser.add_argument("--cell", metavar = 3000, default=3000,
        nargs='?', action = "store", type=int,
        help="Cell count")
    parser.add_argument("--cmo", metavar="cmo.csv",
        nargs='?', action = "store", type=str,
        help="Path to cmo file if applicable")
    parser.add_argument("--feature", metavar="feature.csv",
        nargs='?', action = "store", type=str,
        help="Path to feature barcode reference file if applicable")
    parser.add_argument("--vdjref", metavar="refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0",
        nargs='?', action = "store", type=str,
        help="Path to vdj reference")
    parser.add_argument("--innerprimer", metavar="enrichment_primers.txt",
        nargs='?', action = "store", type=str,
        help="Text file containing one primer per line if non-10x inner enrichment primers were used")
    #parser.add_argument("-f", "--force", action="store_true",
    #    help="Use force-cells flag instead of expect-cells")
    group_cell_number = parser.add_mutually_exclusive_group()
    group_cell_number.add_argument("--force", action="store_true",
        help="Use force-cells flag instead of expect-cells ")
    group_cell_number.add_argument("--expect", action="store_true",
        help="Run Cell Ranger with --expect-cells")
    #parser.add_argument("-i", "--introns", action="store_true",
    #    help="Use include-introns flag")
    ## Since cellranger v7.0.0+, include-introns is true by default
    parser.add_argument("-i", "--exclude_introns", action="store_true",
        help="Use include-introns,false option")
    parser.add_argument("--create_bam", action="store_true",
        help="Use include-introns,false option")
    ## Fix RNA profiling parameters
    parser.add_argument("--probe_set", type=str,
        help="Text files containing probe sets")
    parser.add_argument("--multiplex", type=str,
        help="Text files containing multiplex information")

    args = parser.parse_args(raw_args)
    print(args)

    with open(args.output, 'w', newline="") as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        spamwriter.writerow(['[gene-expression]'])
        spamwriter.writerow(['reference', args.ref])
        ## Parameter for fixed RNA profiling 
        if args.probe_set != None:
            spamwriter.writerow(['probe-set', args.probe_set])
        if args.force:
            spamwriter.writerow(['force-cells', args.cell])
        elif args.expect:
            spamwriter.writerow(['expect-cells', args.cell])
        else:  
            pass
        if args.create_bam:
            spamwriter.writerow(['create-bam', "true"])
        if args.cmo != None:
            spamwriter.writerow(['cmo-set', args.cmo])
        if args.exclude_introns:
            spamwriter.writerow(['include-introns', 'false'])

        if args.feature != None:
            spamwriter.writerow([])
            spamwriter.writerow(['[feature]'])
            spamwriter.writerow(['reference', args.feature])

        if args.vdjref != None:
            spamwriter.writerow([])
            spamwriter.writerow(['[vdj]'])
            spamwriter.writerow(['reference', args.vdjref])
            if args.innerprimer != None:
                spamwriter.writerow(['inner-enrichment-primers', args.innerprimer])

        spamwriter.writerow([])
        multi_crispr = False
        crispr_ref = ""
        with open(args.lib, 'r') as lib:
            line = next(lib)
            for line in lib:
                if "CRISPR Guide Capture" in line:
                    multi_crispr = True
                    line = line.strip().split(',')
                    crispr_ref = line[-1]
            if multi_crispr:
                spamwriter.writerow(["[feature]"])
                spamwriter.writerow(['ref', crispr_ref])
                
            spamwriter.writerow([])
            spamwriter.writerow(['[libraries]'])
            spamwriter.writerow(['fastq_id', 'fastqs', 'lanes', 'feature_types'])
            lib.seek(0)  ## move the file pointer to the begining of the file
            line = next(lib)
            for line in lib:
                line = line.strip().split(',')
                spamwriter.writerow([line[1], line[0], 'Any', line[2]])
        

        if args.cmo != None:
            spamwriter.writerow([])
            spamwriter.writerow(['[samples]'])
            spamwriter.writerow(['sample_id', 'cmo_ids', 'description'])
            with open(args.cmo, 'r') as lib:
                line = next(lib)
                index = 1
                for line in lib:
                    line = line.strip().split(',')
                    spamwriter.writerow(['HTO_%s' % index, line[0]])
                    index += 1
        ## For fixed RNA profiling
        if args.multiplex != None:
            spamwriter.writerow([])
            spamwriter.writerow(['[samples]'])
            spamwriter.writerow(['sample_id', 'probe_barcode_ids', 'description'])
            with open(args.multiplex, 'r') as lib:
                line = next(lib)
                index = 1
                for line in lib:
                    ele = line.strip().split(',')
                    sample_name = re.sub(".csv$", "", args.output)
                    if ele[0] == sample_name:
                        spamwriter.writerow(ele[1:])

if __name__ == '__main__':
    main()
