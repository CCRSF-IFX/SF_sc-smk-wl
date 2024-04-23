import glob, sys
import re
import os

def main(path):
    files = glob.glob('*/outs/molecule_info.h5')
    aggrFile = open('AggregatedDatasets.csv', 'w')

    #aggrFile = open('/Users/chenv3/Downloads/aggregate_samples.csv', 'w')
    # https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/running-pipelines/aggr-aggr
    aggrFile.write('library_id,molecule_h5,cloupe_file,spatial_folder\n')
    files.sort()
    for i in files:
        cloupe_file = os.path.join(path, re.sub("outs/molecule_info.h5", "outs/cloupe.cloupe", i)) 
        spatial_folder = os.path.join(path, re.sub("outs/molecule_info.h5", "outs/spatial", i))
        print([cloupe_file, spatial_folder])
        aggrFile.write('%s,%s/%s,%s,%s\n' % (i.split('/')[0], path, i, cloupe_file, spatial_folder))
    aggrFile.close()

if __name__ == "__main__":
    main(sys.argv[1])
