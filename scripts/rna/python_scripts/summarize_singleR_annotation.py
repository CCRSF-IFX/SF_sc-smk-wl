import os
import sys
import csv
from collections import Counter, defaultdict
import pandas as pd

def find_annotation_files(base_path):
    annotation_files = []
    for root, dirs, files in os.walk(base_path):
        for file in files:
            if file == "annotations_compiled_barcode.csv":
                annotation_files.append(os.path.join(root, file))
    return annotation_files

def summarize_annotations(annotation_files, output_excel):
    # reference -> celltype -> sample -> count
    data = defaultdict(lambda: defaultdict(dict))
    samples = set()
    for file_path in annotation_files:
        sample = os.path.basename(os.path.dirname(os.path.dirname(file_path)))
        samples.add(sample)
        with open(file_path, newline='') as csvfile:
            reader = csv.reader(csvfile)
            header = next(reader)
            references = header[1:]
            counts = {ref: Counter() for ref in references}
            for row in reader:
                for i, ref in enumerate(references, start=1):
                    celltype = row[i]
                    counts[ref][celltype] += 1
            for ref in references:
                for celltype, count in counts[ref].items():
                    data[ref][celltype][sample] = count
    samples = sorted(samples)
    # Write to Excel with one sheet per reference
    with pd.ExcelWriter(output_excel) as writer:
        for ref, celltype_dict in data.items():
            df = pd.DataFrame.from_dict(celltype_dict, orient='index')
            df = df.reindex(columns=samples, fill_value=0)
            df.index.name = 'Cell type'
            df.to_excel(writer, sheet_name=ref)
    print(f"Summary written to {output_excel}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python summarize_singleR_annotation.py <path>")
        sys.exit(1)
    base_path = sys.argv[1]
    annotation_files = find_annotation_files(base_path)
    if not annotation_files:
        print("No annotations_compiled_barcode.csv files found.")
        sys.exit(0)
    output_excel = 'singleR_annotation_summary.xlsx'
    summarize_annotations(annotation_files, output_excel)

if __name__ == "__main__":
    main()
