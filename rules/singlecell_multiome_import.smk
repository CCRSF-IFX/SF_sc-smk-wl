f = open(config.libraries)
line = next(f)
samples = [i.split(',')[0] for i in f.read().strip().split('\n')]
samples = sorted(list(set(samples)))
#print(samples)

f.seek(0)
line = next(f)
fcells = set([i.split(',')[1] for i in f.read().strip().split('\n')])
if len(set(config.unaligned) - set([i for i in config.unaligned for j in fcells if j in i])) > 0:
    raise Exception(f"\nNot all paths in config.unaligned were used in the libraries.csv file. Missing path(s): \n{chr(10).join(set(config.unaligned) - set([i for i in config.unaligned for j in fcells if j in i]))}\n")

f.close()

current_cellranger = program.cellranger_arc
