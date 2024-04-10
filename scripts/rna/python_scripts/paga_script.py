# -*- coding: utf-8 -*-

import os, sys, re, argparse
import pickle
import scanpy as sc
import anndata2ri
from rpy2.robjects import r


def runPAGA(inputFile, h5 = False, workdir="./"):
    if h5:
        adata = sc.read_10x_h5(inputFile)
    else:
        anndata2ri.activate()
        r('library(Seurat)')
        adata = r('as.SingleCellExperiment(readRDS("%s"))' % inputFile)

    os.chdir(workdir)

    adata.var_names_make_unique()

    #Since this is not processed, this part is for clustering
    sc.pp.neighbors(adata)#, n_neighbors=7, n_pcs=20)
    sc.tl.louvain(adata)

    #Running and plotting PAGA
    sc.tl.paga(adata, groups='louvain')
    sc.pl.paga(adata, color=['louvain'], edge_width_scale=0.2, threshold=0.2, save='_louvain')

    #Running and plotting UMAP for comparison
    if 'X_umap' not in adata.obsm_keys():
        sc.tl.umap(adata)
    sc.pl.umap(adata, color=['louvain'], legend_loc='on data', save='_louvain')


    for group in [s for s in adata.obs_keys() if 'snn_res' in s]:
        sc.tl.paga(adata, groups=group)
        #if '%s_colors' % group not in adata.uns_keys():
        #    adata.uns['%s_colors' % group] = sc.pl.palettes.default_102[:len(set(adata.obs[group]))]
        try:
            sc.pl.paga(adata, color=[group], edge_width_scale=0.2, threshold=0.2, save='_%s' % group)
        except:
            print("Couldn't plot paga for %s" % group)
        sc.pl.umap(adata, color=[group], legend_loc='on data', save='_%s' % group)

    if h5:
        f = open('paga.pickle', 'wb')
    else:
        f = open('paga_seurat.pickle', 'wb')
    pickle.dump(adata, f)
    f.close()

def main(raw_args = None):
    parser = argparse.ArgumentParser(description="""Run basic trajectory given either a Seurat RDS object or an h5 file""")
    parser.add_argument("input", metavar="input_file",
                        action="store", type=str,
                        help="Input file: either Seurat RDS or h5 file")
    parser.add_argument("--h5", dest="h5", action="store_true",
                        help="File provided is h5 file")
    parser.add_argument("-w", "--workdir", dest="workdir",
                        action="store", type=str, default="./",
                        help="Working directory to save files in")

    args = parser.parse_args()

    runPAGA(args.input, args.h5, args.workdir)

if __name__=='__main__':

    main()
