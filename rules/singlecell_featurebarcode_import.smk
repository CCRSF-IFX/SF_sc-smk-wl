import re
f = open(config.libraries)
line = next(f)
samples = [i.split(',')[0] for i in f.read().strip().split('\n')]
samples = sorted(list(set(samples)))
#print(samples)

f.seek(0)
line = next(f)
fcells = set([i.split(',')[1] for i in f.read().strip().split('\n')])
#print(fcells)
if len(set(config.unaligned) - set([i for i in config.unaligned for j in fcells if j in i])) > 0:
    raise Exception(f"\nNot all paths in config.unaligned were used in the libraries.csv file. Missing path(s): \n{chr(10).join(set(config.unaligned) - set([i for i in config.unaligned for j in fcells if j in i]))}\n")

f.close()

numcell = getattr(config, "numcells", False) 
dict2 = {}
if numcell != False:
    numcells = numcell.split(',')
    dict2 = dict(zip(samples, numcells))
else:
    dict2 = { item:0 for item in samples} 

def get_cell_number_from_lib_csv():
    flags = ['force-cells', 'expect-cells']
    if hasattr(config, 'libraries'):
        lib_df = pd.read_csv(config.libraries)
        for flag in flags: 
            if flag in lib_df.columns:
                setattr(config, re.sub("-", "", flag), True)
                for index, row in lib_df.iterrows():
                    dict2[row['Name']] = row[flag]

get_cell_number_from_lib_csv()

include_introns = getattr(config, "include_introns", True)

current_cellranger = program.cellranger
