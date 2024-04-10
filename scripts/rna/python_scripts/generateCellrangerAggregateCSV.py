import os, glob, sys

def main(arg1='./', arg2="config_aggr.csv"):
    createAggregateCSV(arg1, arg2)

def createAggregateCSV(arg1, arg2):
    files = glob.glob('*/outs/molecule_info.h5')
    files.sort()
    path = os.getcwd()
    f = open(arg1 + '/' + arg2, 'w')
    f.write("library_id,molecule_h5\n")
    for i in files:
        f.write(','.join([i.split('/')[0], path+'/'+i]))
        f.write('\n')

if __name__ == "__main__":
    if len(sys.argv) == 1:
        main()
    elif len(sys.argv) == 2:
        main(sys.argv[1])
    elif len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        print("Expecting only one or two input variables, taking first two variables for processing.")
        main(sys.argv[1], sys.argv[2])
