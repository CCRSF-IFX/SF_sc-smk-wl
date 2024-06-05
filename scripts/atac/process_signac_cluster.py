import argparse


def parse_file(file_name, out_dir):
    f = open(file_name)
    clusters = dict()
    for line in f:
        line = line.strip().split(',')
        barcodes = clusters.get(line[1], list())
        barcodes.append(line[0])
        clusters[line[1]] = barcodes
    f.close()

    for i in clusters:
        f = open("%s/Cluster%s_cells.txt" % (out_dir, i), "w")
        for cell in clusters[i]:
            f.write("CB:Z:%s\n" % cell)
        f.close()


def main(raw_args = None):
    parser = argparse.ArgumentParser(description="Parse signac cluster file for bam file generation")
    #parser.add_argument('sample_name', metavar="Sample", action="store",
    #    type=str, help="Name of sample")
    parser.add_argument('file_name', metavar="clusters.csv", action="store",
        type=str, help="CSV file containing signac cluster results")
    parser.add_argument('out_dir', metavar="output_directory", action="store",
        type=str, nargs='?', const=1, default="./", help="Output directory for cluster results")

    args = parser.parse_args(raw_args)
    parse_file(args.file_name, args.out_dir)

if __name__ == '__main__':
    main()
