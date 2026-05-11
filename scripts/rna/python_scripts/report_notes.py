def get_wordreport_notes(application):
    if "RNA" not in application.upper():
        return None

    return (
        "Number of Cell - The number of barcodes associated with cell-containing partitions.\n"
        "Mean Reads per Cell - The total number of sequenced reads divided by the estimated number of cells.\n"
        "Median Genes per Cell - The median number of genes detected per cell barcode.\n"
        "Total Genes Detected - The number of genes detected.\n"
        "Median UMI Counts per Cell - The median number of UMI counts per cell barcode."
    )
