#!/usr/bin/env Rscript

option_list <- list(
  optparse::make_option(c("--rds"),
        type = "character", default = NULL,
        help = "A list of RDS files of Seurat objects", metavar = "character"
  ),
  optparse::make_option(c("--sampleid"),
    type = "character", default = "scRNA",
    help = "Sample ID", metavar = "character"
  ),
  optparse::make_option(c("--genome"),
    type = "character", default = "hg38",
    help = "Genome build [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--markerList"),
    type = "character", default = NULL,
    help = "Marker gene list", metavar = "character"
  ),
  optparse::make_option(c("--outdir"),
    type = "character", default = "output/",
    help = "Output RDS file name [default= %default]", metavar = "character"
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)
system(paste0("mkdir -p ", opt$outdir))




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
  tryCatch({
    p <- plotDeltaDistribution(pred, show="delta.med", ncol=3)
    print(p)
  }, error=function(e) message("Skipping delta distribution plot: ", conditionMessage(e)))
  dev.off()
  seur@meta.data[name] <- ''
  seur@meta.data[rownames(pred),name] <- pred$pruned.labels
  Idents(seur) <- seur@meta.data[,name]
  png(paste0("UMAP_", name, ".png"))
  print(DimPlot(seur))
  dev.off()
  png(paste0("TSNE_", name, ".png"))
  print(DimPlot(seur, reduction='tsne'))
  dev.off()
  dir.create("gene_list_plots")
  for (i in 1:dim(markers)[1]) {
    genes.use <- as.character(markers[row.names(markers)[i],])
    numGenes <- length(genes.use[genes.use %in% row.names(seur)])
    if (sum(genes.use %in% row.names(seur)) > 1) {
      png(paste0('gene_list_plots/ViolinPlot_', name, '_', row.names(geneList)[i],'.png'), height = 5, width = (length(genes.use[genes.use %in% row.names(seur)])+1)*5, units='in', res=300)
      print(suppressWarnings(VlnPlot(seur, genes.use[genes.use %in% row.names(seur)], ncol=numGenes)))
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
      print(suppressWarnings(FeaturePlot(seur, genes.use[genes.use %in% row.names(seur)], ncol=numGenes)))
      dev.off()
    }
  }
}

# The objects returned by functions:  ImmGenData, MouseRNAseqData, 
# HumanPrimaryCellAtlasData and BlueprintEncodeData are saved under 
# /home/docker/. So you only need to load the RDS files directly. 
if (startsWith(opt$genome, "mm10") | startsWith(opt$genome, "mm39")) {
  #immgen <- ImmGenData()# - mouse
  #mouserna <- MouseRNAseqData()# - mouse
  immgen <- readRDS("/home/docker/ImmGenData.rds")
  mouserna <- readRDS("/home/docker/MouseRNAseqData.rds")
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
} else if (startsWith(opt$genome, "hg38") | startsWith(opt$genome, "hg19")) {
  #hpca <- HumanPrimaryCellAtlasData()
  #blueprint <- BlueprintEncodeData()
  hpca <- readRDS("/home/docker/HumanPrimaryCellAtlasData.rds")
  blueprint <- readRDS("/home/docker/BlueprintEncodeData.rds")
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

suppressWarnings(generateMarkerPlots(seur, geneList))
saveRDS(seur, 'seur_10x_cluster_singler.rds')

sessionInfo()



