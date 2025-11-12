#!/usr/bin/env Rscript

# Adapted from cellsnake here: https://github.com/sinanugur/scrna-workflow

option_list <- list(
  optparse::make_option(c("--min.cells"),
    type = "integer", default = 1,
    help = "Min cells [default= %default]", metavar = "integer"
  ),
  optparse::make_option(c("--min.features"),
    type = "integer", default = 1,
    help = "Min features, nFeature_RNA [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--data.dir"),
    type = "character", default = NULL,
    help = "Data directory", metavar = "character"
  ),
  optparse::make_option(c("--sampleid"),
    type = "character", default = "scRNA",
    help = "Sample ID", metavar = "character"
  ),
  optparse::make_option(c("--genome"),
    type = "character", default = "hg38",
    help = "Genome build [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--outdir"),
    type = "character", default = "output/",
    help = "Output RDS file name [default= %default]", metavar = "character"
  )
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

system(paste0("mkdir -p ", opt$outdir))




library(Seurat)
library(Matrix)
library(MASS)
library(dplyr)
library(reshape2)
library(ggplot2)
library(URD)
library(cluster)


opt                                                              

system(paste0("mkdir -p ", opt$outdir))
setwd(opt$outdir) 

# Initialize variable (optional but explicit)
count_mtx <- NULL  

# Define file paths
parsebio_files <- c("count_matrix.mtx", "all_genes.csv", "cell_metadata.csv")
tenx_files <- c("matrix.mtx", "matrix.mtx.gz", 
                "features.tsv", "features.tsv.gz",
                "genes.tsv", "genes.tsv.gz",
                "barcodes.tsv", "barcodes.tsv.gz")

# Check file existence
has_parsebio <- all(file.exists(file.path(opt$data.dir, parsebio_files)))
has_10x <- any(file.exists(file.path(opt$data.dir, tenx_files)))

# Select reader
if (has_parsebio) {
  message("Detected ParseBio format.")
  count_mtx <- ReadParseBio(data.dir = opt$data.dir)
} else if (has_10x) {
  message("Detected 10x Genomics format.")
  count_mtx <- Read10X(data.dir = opt$data.dir)
} else {
  stop("No recognized data format found in: ", opt$data.dir)
}


seur <- CreateSeuratObject(counts = count_mtx,
                                 min.cells = opt$min.cells,
                                 min.features = opt$min.features,
                                 project = opt$sampleid)

rm(count_mtx)

#find mitochondrial genes
if (startsWith(opt$genome, "mm10") | startsWith(opt$genome, "mm39")) {
        seur[["percent.mito"]] <- PercentageFeatureSet(seur, pattern="^mt-")
} else if (startsWith(opt$genome, "hg19") | startsWith(opt$genome, "hg38")) {
        seur[["percent.mito"]] <- PercentageFeatureSet(seur, pattern="^MT-")
}

head(seur@meta.data)


png("VlnPlot_PreFilter.png", height=7, width=7, units='in', res=200)
VlnPlot(seur, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

plot1 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png("FeatureScatter_PreFilter.png", height=7, width=10, units='in', res=200)
plot1 | plot2
dev.off()

cS <- data.frame(libSize=seur$nCount_RNA, geneDetect=seur$nFeature_RNA)
p_hi <- 1e-3 #p-value for filtering doublets
p_lo <- 1e-2 #p-value for filtering poor libraries
fitLibSize <- fitdistr(cS$libSize,"negative binomial")
umi.upper.limit <- qnbinom(p_hi,size=fitLibSize$estimate["size"],
                 mu=fitLibSize$estimate["mu"],lower.tail=F)
umi.lower.limit <- qnbinom(p_lo,size=fitLibSize$estimate["size"],
                 mu=fitLibSize$estimate["mu"],lower.tail=T)
fitGeneDetect <- fitdistr(cS$geneDetect,"negative binomial")
gene.upper.limit <- qnbinom(p_hi,size=fitGeneDetect$estimate["size"],
                 mu=fitGeneDetect$estimate["mu"],lower.tail=F)
gene.lower.limit <- qnbinom(p_lo,size=fitGeneDetect$estimate["size"],
                 mu=fitGeneDetect$estimate["mu"],lower.tail=T)

cur.mad <- mad(seur$percent.mito)
cur.med <- median(seur$percent.mito)
diff.val <- 4 * cur.mad
mito.upper.limit <- cur.med + diff.val

temp_doublets <- (cS$libSize > umi.upper.limit) | (cS$geneDetect > gene.upper.limit) #doublets IDed based on high library size or genes detected
temp_crapLibs <- (cS$libSize < umi.lower.limit) | (cS$geneDetect < gene.lower.limit) #poor libraries IDed based on low library size or genes detected

print(mito.upper.limit)

#filter_summary <- t(data.frame(c("Doublets"=sum(temp_doublets), "Poor-Quality"=sum(temp_crapLibs), "Mitochondrial"=sum(seur$percent.mito > mito.upper.limit), "Total Filtered"=sum(temp_doublets | temp_crapLibs | seur$percent.mito > mito.upper.limit))))

filter_summary <- t(data.frame(c(
  "Doublets"=sum(temp_doublets),
  "Poor-Quality"=sum(temp_crapLibs),
  "Mitochondrial"=sum(seur$percent.mito >= mito.upper.limit),
  "Total Filtered"=sum(temp_doublets | temp_crapLibs | seur$percent.mito >= mito.upper.limit)
)))

dim(seur)

rownames(filter_summary) = opt$sampleid
filter_summary

write.table(filter_summary, 'FilterNumbers.csv', sep=',', quote=FALSE, row.names=FALSE)

seur <- subset(seur, subset = nFeature_RNA >= gene.lower.limit & nFeature_RNA <= gene.upper.limit)
seur <- subset(seur, subset = nCount_RNA   >= umi.lower.limit  & nCount_RNA   <= umi.upper.limit)
# Only apply mitochondrial filtering if the threshold is positive
if (mito.upper.limit > 0) {
  seur <- subset(seur, subset = percent.mito <= mito.upper.limit)
}

# seur@misc$qc_cutoffs <- list(
#   gene = c(lower = gene.lower.limit, upper = gene.upper.limit),
#   umi  = c(lower = umi.lower.limit,  upper = umi.upper.limit),
#   mito = mito.upper.limit
# )

write.table(data.frame(umi.upper.limit, umi.lower.limit, gene.upper.limit, gene.lower.limit, mito.upper.limit), 
            'FilterThresholds.csv', sep=',', quote=FALSE, row.names = FALSE)

png("VlnPlot_Filtered.png", height=7, width=7, units='in', res=200)
VlnPlot(seur, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

#new plot
plot1 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png("FeatureScatter_PostFilter.png", height=7, width=10, units='in', res=300)
plot1 | plot2
dev.off()

seur <- SCTransform(seur, vars.to.regress = "percent.mito", return.only.var.genes = FALSE, verbose = FALSE)

saveRDS(seur, file = "seur_10x_preprocessed_object.rds")

###URD
if (dim(seur)[[2]] > 10000) {
  numPCs <- 50
} else if (dim(seur)[[2]] < 500) {
  numPCs <- 10
} else {
  #inputTags <- as.matrix(read.csv(expressionFile, row.names = 1))
  mat1 <- as.matrix(GetAssayData(seur))
  cat("createURD...")
  test <- suppressWarnings(createURD(count.data = mat1, min.cells=3, min.counts=3)) # )
  cat("calcPCA...")
  test <- suppressWarnings(calcPCA(test, mp.factor = 2))
  write.table(test@pca.sig,"URD.txt")
  png("URD.png", height=7, width=7, units='in', res=300)
  pcSDPlot(test)
  dev.off()

  numPCs <- max(10, sum(test@pca.sig))
}

top10 <- head(VariableFeatures(seur), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seur)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
png('VariableFeatures.png', height=7, width=10, units='in', res=300)
plot1 | plot2
dev.off()

seur <- RunPCA(seur, 
               features = VariableFeatures(seur), 
               npcs=max(20,numPCs), do.print = TRUE, 
               ndims.print = 1:5, nfeatures.print = 5)

pdf("VizPCAPlot.pdf")
for (i in  seq(1,20,2)) {
  j = i + 1
  print(VizDimLoadings(seur, i:j))
}
dev.off()

#plot PCA
pdf("AllPCAPlot.pdf")
for (i in c(1:10)) {
  print(DimPlot(seur, dims=c(i, i+1), reduction="pca"))
}
dev.off()

#PC heatmap
pdf("PC_HeatmapPlot.pdf")
for (i in c(1:10)) {
  DimHeatmap(seur, dims = i, cells = 500, balanced = TRUE)
}
dev.off()

#make PC elbow plot
pdf("PC_ElbowPlot.pdf")
ElbowPlot(seur)
dev.off()

resolutions <- c(0.1, 0.3, 0.6, 0.8)

seur <- FindNeighbors(seur, dims=1:numPCs)
seur <- RunTSNE(seur, dims=1:numPCs)
write.csv(Embeddings(seur, reduction='tsne'), file = "tSNECoordinates.csv")
seur <- RunUMAP(seur, dims=1:numPCs)
write.csv(Embeddings(seur, reduction='umap'), file = "UMAPCoordinates.csv")

runRes <- c()
tsnePlots <- list()
umapPlots <- list()
for (res in resolutions) {
  seur <- FindClusters(seur, dims=1:numPCs, resolution = res, print.output = 0, save.SNN = T)
  tsne <- TSNEPlot(seur) + ggtitle(paste(numPCs,"PCs_res", res, sep="")) +
    theme(plot.title = element_text(hjust = 0.5))
  tsnePlots[[as.character(res)]] <- tsne
  png(paste("TSNEPlotwith",numPCs,"PCs_", res, ".png", sep=""), height=7, width=7, units='in', res=300)
  print(tsne)
  dev.off()
  umap <- DimPlot(seur, reduction="umap") + ggtitle(paste(numPCs,"PCs_res", res, sep="")) +
    theme(plot.title = element_text(hjust = 0.5))
  umapPlots[[as.character(res)]] <- umap
  png(paste("UMAPPlotwith",numPCs,"PCs_", res, ".png", sep=""), height=7, width=7, units='in', res=300)
  print(umap)
  dev.off()

  try({
    seur.markers <- FindAllMarkers(object = seur, logfc.threshold = 0.25, only.pos=TRUE)
    write.csv(seur.markers %>% group_by(cluster) %>% top_n(-100,
                                                           p_val), paste("top100markers_pc", numPCs, "_res", res, ".csv", sep = ""))
    saveRDS(seur.markers, paste("markers_res", res, ".rds", sep = ""))
    runRes <- append(runRes, res)})
}

#save object
saveRDS(seur, file = "seur_10x_cluster_object.rds")

pdf("TSNEPlots.pdf")
for (res in tsnePlots){
  print(res)
}
dev.off()

pdf("UMAPPlots.pdf")
for (res in umapPlots){
  print(res)
}
dev.off()

##Create Silhoutte Plots
for (res in runRes){
  coord <- Embeddings(seur, reduction.type='pca')[,1:numPCs]
  Idents(seur) <- seur@meta.data[[paste0('SCT_snn_res.', res)]]
  clusters <- Idents(seur)
  d <- dist(coord, method="euclidean")
  sil<-silhouette(as.numeric(clusters), dist=d)
  #silPlot <- recordPlot()
  pdf(paste0("SilhouettePlot_res",res,".pdf"))#, height=7, width=7, units='in', res=300)
  plot(sil, col=as.factor(clusters[order(clusters, decreasing=FALSE)]), main=paste("Silhouette plot of Seurat clustering - resolution ", res, sep=""), lty=2)
  abline(v=mean(sil[,3]), col="red4", lty=2)
  dev.off()
}

##Remove resolutions that failed marker generation
print(runRes)
for (res in setdiff(resolutions, runRes)){
  seur@meta.data[paste0('SCT_snn_res.', res)] <- NULL
}

write(min(runRes), "minRes.txt")

sessionInfo()



