import os, argparse
import pickle
import scvelo as scv
import scanpy as sc
import anndata2ri
import scmomentum as scm
from rpy2.robjects import r
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def runSCVelo(loom_input, seurat_file, workdir="./", sampleName = None):
    adata_full = sc.read(loom_input)
    adata_full.var_names_make_unique()

    anndata2ri.activate()
    r('library(Seurat)')
    seurat = r('as.SingleCellExperiment(readRDS("%s"))' % seurat_file)

    if sampleName == None:
        sampleName = seurat.obs['orig.ident'][0]

    os.chdir(workdir)

    names = [i.split(':')[1][:-1] for i in adata_full.obs_names]
    adata = adata_full[[sampleName + ':%sx' % name for name in seurat.obs_names if name in names] ,[name for name in seurat.var_names if name in adata_full.var_names]]

    adata.obs['SCT_snn_res.0.6'] = seurat.obs['SCT_snn_res.0.6'].add_prefix('%s:' % sampleName).add_suffix('x')
    adata.obsm = seurat.obsm

    scv.pp.filter_genes(adata, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata)
    scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
    scv.pp.log1p(adata)

    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode='dynamical')
    scv.tl.velocity_graph(adata)

    scv.pl.velocity_embedding_stream(adata, basis='umap', color='SCT_snn_res.0.6', save='stream')
    scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=2, dpi=120, color='SCT_snn_res.0.6', save='vector')

    scv.tl.velocity_confidence(adata)
    keys = 'velocity_length', 'velocity_confidence'
    scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save='length_confidence')

    df = adata.var
    df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]

    fig2 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=3, nrows=1, figure=fig2)
    f2_ax1 = fig2.add_subplot(spec2[0, 0])
    f2_ax2 = fig2.add_subplot(spec2[0, 1])
    f2_ax3 = fig2.add_subplot(spec2[0, 2])
    f2_ax1.hist(df['fit_alpha'])
    f2_ax1.set_xlabel('transcription rate')
    f2_ax2.hist(df['fit_beta'] * df['fit_scaling'])
    f2_ax2.set_xlabel('splicing rate')
    f2_ax3.hist(df['fit_gamma'])
    f2_ax3.set_xlabel('degradation rate')
    fig2.tight_layout()
    fig2.set_size_inches(18, 6)
    fig2.savefig('kinetic_rate.png', dpi=128, bbox_inches='tight', pad_inches=0.0)

    scv.tl.latent_time(adata)
    scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save='latent_time')

    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
    scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='SCT_snn_res.0.6', n_convolve=100, save='latent_time_top_genes')

    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
    scv.pl.scatter(adata, basis=top_genes[:15], color='SCT_snn_res.0.6', ncols=5, frameon=False, save='top_likelihood_genes')

    scv.tl.rank_dynamical_genes(adata, groupby='SCT_snn_res.0.6')
    df = scv.get_df(adata, 'rank_dynamical_genes/names')
    df.head(5)

    for cluster in  df.columns:
        scv.pl.scatter(adata, df[cluster][:5], ylabel=cluster, frameon=False, color='SCT_snn_res.0.6', save='top_likelihood_genes_cluster' + cluster)

    f = open('scvelo_dynamic.pickle', 'wb')
    pickle.dump(adata, f)
    f.close()


    group = 'SCT_snn_res.0.6'
    adata = scm.network_inference.preprocess(adata,group)

    clusters = list(set(adata.obs[group]))
    mode = 'highexp'
    size = 100

    for cluster in clusters:
        scm.network_inference.predict_network(adata,cluster,mode,size,group)

    for cluster in clusters:
        adata.uns[mode+'-'+str(size)][cluster].to_csv('scMomentum_highexp_100_cluster%s.csv' % cluster)

    f = open('scmomentum.pickle', 'wb')
    pickle.dump(adata, f)
    f.close()

def main(raw_args = None):
    parser = argparse.ArgumentParser(description="""Run scVelo and scMomentum given a velocyto loom file and Seurat RDS object""")
    parser.add_argument("loom", metavar="loom_file",
                        action="store", type=str,
                        help="Input velocyto loom file")
    parser.add_argument("seurat_rds", metavar="seurat_file",
                        action="store", type=str,
                        help="Input Seurat RDS file")
    parser.add_argument("-w", "--workdir", dest="workdir",
                        action="store", type=str, default="./",
                        help="Working directory to save files in")
    parser.add_argument("-n", dest="sample_name",
                        action="store", type=str, default=None,
                        help="Sample name")

    args = parser.parse_args()

    print(args)
    runSCVelo(args.loom, args.seurat_rds, args.workdir, args.sample_name)

if __name__=='__main__':

    main()
