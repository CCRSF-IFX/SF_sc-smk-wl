script_path = getwd()                                # %exclude_jupyterlab%

script_path                                          # %exclude_jupyterlab%

system(paste0('Rscript ', script_path,                # %exclude_jupyterlab%
        '/sc_singleR_opt.R --genome="hg38" --markerList="/Volumes/ccrsf-ifx/Software/scripts/bin/currentsnake/single_cell/gene_lists/human_gene_list.csv" --outdir="/Volumes/ccrsf-static/Analysis/xies4/github_repos/pipeline_dev_test/singleR" --rds="/Volumes/ccrsf-static/Analysis/xies4/github_repos/pipeline_dev_test/test_dir/seur_10x_cluster_object.rds"'),  # %exclude_jupyterlab%
       intern = T) # %exclude_jupyterlab%

opt = readRDS("/Volumes/ccrsf-static/Analysis/xies4/github_repos/pipeline_dev_test/singleR/opt.rds")         # %exclude_jupyterlab%

opt

system(paste0("mkdir -p ", opt$outdir))
setwd(opt$outdir) 

library(Seurat)
library(SingleR)
library(celldex)

seur <- readRDS(opt$rds)

sce <- as.SingleCellExperiment(seur)


sce <- scater::logNormCounts(sce)

geneList <- read.csv(opt$markerList, header=FALSE, row.names=1, stringsAsFactors = FALSE)

## plotScoreDistribution: use this plot to check that outlier detection in pruneScores() behaved sensibly. 
## Labels with especially low deltas may warrant some additional caution in their interpretation.
generatePlots <- function(name, pred, seur, markers){
  saveRDS(pred, paste0('pred_', name, '.rds'))
  #pdf(paste0('heatmap_', name, '.pdf'))
  #print(plotScoreHeatmap(pred, clusters=sce$res.0.1))
  #dev.off()
  pdf(paste0('ScoreDistribution_', name, '.pdf'), height=15)
  # plotScoreDistribution deprecated
  #print(plotScoreDistribution(pred, show = "delta.med", ncol = 3, show.nmads = 3))
  print(plotDeltaDistribution(pred, show = "delta.med", ncol = 3))
  dev.off()
  seur@meta.data[name] <- ''
  seur@meta.data[rownames(pred),name] <- pred$pruned.labels
  Idents(seur) <- seur@meta.data[,name]
  pdf(paste0("UMAP_", name, ".pdf"))
  print(DimPlot(seur))
  dev.off()
  pdf(paste0("TSNE_", name, ".pdf"))
  print(DimPlot(seur, reduction='tsne'))
  dev.off()
  dir.create("gene_list_plots")
  for (i in 1:dim(markers)[1]) {
    genes.use <- as.character(markers[row.names(markers)[i],])
    numGenes <- length(genes.use[genes.use %in% row.names(seur)])
    if (sum(genes.use %in% row.names(seur)) > 1) {
      png(paste0('gene_list_plots/ViolinPlot_', name, '_', row.names(geneList)[i],'.png'), height = 5, width = (length(genes.use[genes.use %in% row.names(seur)])+1)*5, units='in', res=300)
      print(VlnPlot(seur, genes.use[genes.use %in% row.names(seur)], ncol=numGenes))
      dev.off()
    }
  }
  return(seur)
}

generateMarkerPlots <- function(seur, markers) {
  dir.create("gene_list_plots")
  for (i in 1:dim(markers)[1]) {
    genes.use <- as.character(markers[row.names(markers)[i],])
    numGenes <- length(genes.use[genes.use %in% row.names(seur)])
    if (sum(genes.use %in% row.names(seur)) > 1) {
      png(paste0('gene_list_plots/FeaturePlot_', row.names(geneList)[i],'.png'), height = 5, width = (length(genes.use[genes.use %in% row.names(seur)])+1)*5, units='in', res=300)
      print(FeaturePlot(seur, genes.use[genes.use %in% row.names(seur)], ncol=numGenes))
      dev.off()
    }
  }
}

# The objects returned by functions:  ImmGenData, MouseRNAseqData, 
# HumanPrimaryCellAtlasData and BlueprintEncodeData are saved under 
# /home/docker/. So you only need to load the RDS files directly. 
if (opt$genome == "mm10") {
  #immgen <- ImmGenData()# - mouse
  #mouserna <- MouseRNAseqData()# - mouse
  immgen <- readRDS("~/ImmGenData.rds")
  mouserna <- readRDS("~/MouseRNAseqData.rds")
  pred.multi <- SingleR(test = sce,
      ref = list(IG=immgen, MRNA=mouserna),
      labels = list(immgen$label.main, mouserna$label.main))
  seur <- generatePlots('immgen_mourserna', pred.multi, seur, geneList)
  temp <- table(pred.multi$pruned.labels)
  write.csv(sort(temp, decreasing=TRUE), 'annotations_Immgen-MouseRNA_general.csv', row.names=FALSE)

  pred.immgen <- SingleR(test = sce,
      ref = immgen, labels = immgen$label.main)
  seur <- generatePlots('immgen', pred.immgen, seur, geneList)
  temp <- table(pred.immgen$pruned.labels)
  write.csv(sort(temp, decreasing=TRUE), 'annotations_Immgen_general.csv', row.names=FALSE)

  pred.mouserna <- SingleR(test = sce,
      ref = mouserna, labels = mouserna$label.main)
  seur <- generatePlots('mourserna', pred.mouserna, seur, geneList)
  temp <- table(pred.mouserna$pruned.labels)
  write.csv(sort(temp, decreasing=TRUE), 'annotations_MouseRNA_general.csv', row.names=FALSE)

  write.csv(seur@meta.data[,c('immgen_mourserna', 'immgen', 'mourserna')], 'annotations_compiled_barcode.csv', quote=FALSE)
} else if (opt$genome == "hg38") {
  #hpca <- HumanPrimaryCellAtlasData()
  #blueprint <- BlueprintEncodeData()
  hpca <- readRDS("~/HumanPrimaryCellAtlasData.rds")
  blueprint <- readRDS("~/BlueprintEncodeData.rds")
  pred.multi <- SingleR(test = sce,
      ref = list(BP=blueprint, HPCA=hpca),
      labels = list(blueprint$label.main, hpca$label.main))
  seur <- generatePlots('blueprintencode_hpca', pred.multi, seur, geneList)
  temp <- table(pred.multi$pruned.labels)
  write.csv(sort(temp, decreasing=TRUE), 'annotations_BlueprintENCODE-HPCA_general.csv', row.names=FALSE)

  pred.blueprint <- SingleR(test = sce,
      ref = blueprint, labels = blueprint$label.main)
  seur <- generatePlots('blueprintencode', pred.blueprint, seur, geneList)
  temp <- table(pred.blueprint$pruned.labels)
  write.csv(sort(temp, decreasing=TRUE), 'annotations_BlueprintENCODE_general.csv', row.names=FALSE)

  pred.hpca <- SingleR(test = sce,
      ref = hpca, labels = hpca$label.main)
  seur <- generatePlots('hpca', pred.hpca, seur, geneList)
  temp <- table(pred.hpca$pruned.labels)
  write.csv(sort(temp, decreasing=TRUE), 'annotations_HPCA_general.csv', row.names=FALSE)

  write.csv(seur@meta.data[,c('blueprintencode_hpca', 'blueprintencode', 'hpca')], 'annotations_compiled_barcode.csv', quote=FALSE)
}

generateMarkerPlots(seur, geneList)
saveRDS(seur, 'seur_10x_cluster_singler.rds')


sessionInfo()

notebook_prefix = "sc_singleR"                                                         # %exclude_jupyterlab%
notebook_name = paste0(notebook_prefix, ".ipynb")                                     # %exclude_jupyterlab%
notebook_r = paste0(script_path, "/", paste0(notebook_prefix, ".r"))                  # %exclude_jupyterlab%
notebook_path = paste0(script_path, "/", notebook_name)                               # %exclude_jupyterlab%
opt_name = paste0(script_path, "/", sub(".ipynb", "_opt.R", notebook_name))           # %exclude_jupyterlab%
output = paste0(script_path, "/", sub(".ipynb", ".prod.R", notebook_name))            # %exclude_jupyterlab%
cmd1 = paste0("jupyter nbconvert --to script --output ",                              # %exclude_jupyterlab%
             notebook_prefix, ' ', notebook_path, "> /dev/null 2>&1 ")                # %exclude_jupyterlab%
cmd1                                                                                  # %exclude_jupyterlab%
system(cmd1, intern = TRUE)                                                            # %exclude_jupyterlab%

cmd2 = paste0('cat ', opt_name, ' ', notebook_r,                                      # %exclude_jupyterlab%
             ' |grep -v exclude_jupyterlab > ', output,  ' 2>&1')                     # %exclude_jupyterlab%
cmd2                                                                                  # %exclude_jupyterlab%
system(cmd2, intern = T)                                                              # %exclude_jupyterlab%
system(paste0("rm ", notebook_r))                                                     # %exclude_jupyterlab%  


