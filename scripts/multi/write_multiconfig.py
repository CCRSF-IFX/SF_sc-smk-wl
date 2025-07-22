#!/usr/bin/env python

import argparse, csv, re, os

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
    parser.add_argument("--hashedabc", action="store_true",
        help="Hashing with Antibody Capture libraries")
    parser.add_argument("--features", metavar="feature.csv",
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
        help="Use create_bam,true option")
    parser.add_argument("--disable_lib_check", action="store_true",
        help="Use check-library-compatibility,false option (default option)")
    ## Fix RNA profiling parameters
    parser.add_argument("--probe_set", type=str,
        help="Text files containing probe sets")
    parser.add_argument("--multiplex", type=str,
        help="Text files containing multiplex information")
    

    args = parser.parse_args(raw_args)
    print(args)
    import pandas as pd
    df = pd.read_csv(args.lib, header=0)
    if 'ocm_barcode_ids' in df.columns:
        if df['ocm_barcode_ids'].any():
            print("OCM platform detected, using OCM specific libraries file format")
            args.ocm = True
        else:
            args.ocm = False

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
        if args.disable_lib_check:
            spamwriter.writerow(['check-library-compatibility', "false"])
        if args.cmo != None and args.hashedabc != True:
            spamwriter.writerow(['cmo-set', args.cmo])
        if args.exclude_introns:
            spamwriter.writerow(['include-introns', 'false'])

        if args.features != None:
            spamwriter.writerow([])
            spamwriter.writerow(['[feature]'])
            spamwriter.writerow(['reference', args.features])

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
            if args.hashedabc != True:
                spamwriter.writerow(['sample_id', 'cmo_ids', 'description'])
                with open(args.cmo, 'r') as lib:
                    line = next(lib)
                    index = 1
                    for line in lib:
                        line = line.strip().split(',')
                        spamwriter.writerow(['HTO_%s' % index, line[0], line[0]])
                        index += 1
            else: ## hashed antibody capture
                with open(args.cmo, 'r') as lib:
                    index = 1
                    for line in lib:
                        line = line.strip().split(',')
                        spamwriter.writerow(line)
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
        if args.ocm:
            samples = ','.join(df['sample_id'].unique()).split(',')
            
            spamwriter.writerow([])
            spamwriter.writerow(['[samples]'])
            header = ['sample_id', 'ocm_barcode_ids', 'description']
            record_cell_number = {}
            record_ocm_ids = {}
            ocm_ids = ','.join(df['ocm_barcode_ids'].unique()).split(',')
            if len(samples) != len(ocm_ids):
                raise ValueError(f"Number of samples ({len(samples)}) does not match number of ocm_ids ({len(ocm_ids)})")
            # Map sample_id to ocm_id
            record_ocm_ids = dict(zip(samples, ocm_ids))

            # Ensure expected_cells and samples have the same number
            if 'expect_cells' in df.columns:
                expected_cells = ','.join(df['expect_cells'].unique()).split(',')
                if len(samples) != len(expected_cells):
                    raise ValueError(f"Number of samples ({len(samples)}) does not match number of expected_cells ({len(expected_cells)})")
                # Map sample_id to cell number
                record_cell_number = dict(zip(samples, expected_cells))
                header.append('expect_cells')
                header = ['sample_id', 'ocm_barcode_ids', 'description', 'expect_cells']
            if "force_cells" in df.columns:
                force_cells = ','.join(df['force_cells'].unique()).split(',')
                if len(samples) != len(force_cells):
                    raise ValueError(f"Number of samples ({len(samples)}) does not match number of force_cells ({len(force_cells)})")
                # Map sample_id to cell number
                record_cell_number = dict(zip(samples, force_cells))
                header.append('force_cells')
                header = ['sample_id', 'ocm_barcode_ids', 'description', 'force_cells']
            spamwriter.writerow(header)
            for sample in samples:
                #A name to identify a multiplexed sample. Must be alphanumeric with hyphens 
                #and/or underscores, and less than 64 characters.
                if len(sample) > 64:
                    raise ValueError(f"Sample name '{sample}' exceeds 64 characters.")
                row = [sample]
                lib_name = os.path.basename(args.lib).replace('_libraries.csv', '')
                print(lib_name)
                ocm_ids = df[df['Sample'] == lib_name]['ocm_barcode_ids'].unique()
                print(df[df['Sample'] == sample])
                row.append(record_ocm_ids[sample])
                row.append(sample)  # description is the sample name
                if 'expect_cells' in df.columns or 'force_cells' in df.columns:
                    row.append(record_cell_number[sample])
                spamwriter.writerow(row)

if __name__ == '__main__':
    main()
