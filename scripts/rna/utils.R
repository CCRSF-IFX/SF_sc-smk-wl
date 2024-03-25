#' Filter Seurat object
#'
#' @param seur Seurat object.
#' @param p_hi Upper end pvalue cutoff for filtering.
#' @param p_hi Lower end pvalue cutoff for filtering.
#' @return A list.
#' @examples

FilterSeuratObject <- function(seur, p_hi = 1e-3, p_lo = 1e-2){
  # p_hi <- 1e-3 #p-value for filtering doublets
  # p_lo <- 1e-2 #p-value for filtering poor libraries
  cS <- data.frame(libSize=seur@meta.data$nCount_RNA, geneDetect=seur@meta.data$nFeature_RNA)

  fitLibSize <- MASS::fitdistr(cS$libSize,"negative binomial")
  umi.upper.limit <- qnbinom(p_hi,size=fitLibSize$estimate["size"],
                             mu=fitLibSize$estimate["mu"],lower.tail=F)
  umi.lower.limit <- qnbinom(p_lo,size=fitLibSize$estimate["size"],
                             mu=fitLibSize$estimate["mu"],lower.tail=T)
  ## only around 0.5% of the cells will be filtered
  perc <- 100*sum(seur@meta.data$nCount_RNA < umi.upper.limit)/length(seur@meta.data$nCount_RNA)
  cat("The percentage of cells will be filterd by umi.upper.limit is: ", perc, "\n")
  fitGeneDetect <- MASS::fitdistr(cS$geneDetect,"negative binomial")
  gene.upper.limit <- qnbinom(p_hi,size=fitGeneDetect$estimate["size"],
                              mu=fitGeneDetect$estimate["mu"],lower.tail=F)
  gene.lower.limit <- qnbinom(p_lo,size=fitGeneDetect$estimate["size"],
                              mu=fitGeneDetect$estimate["mu"],lower.tail=T)
  temp_doublets <- (cS$libSize > umi.upper.limit) | (cS$geneDetect > gene.upper.limit) #doublets IDed based on high library size or genes detected
  temp_crapLibs <- (cS$libSize < umi.lower.limit) | (cS$geneDetect < gene.lower.limit) #poor libraries IDed based on low library size or genes detected
  gene.upper.limit
  gene.lower.limit

  # median absolute deviation
  cur.mad <- mad(seur@meta.data$percent.mt)
  cur.med <- median(seur@meta.data$percent.mt)
  diff.val <- 4 * cur.mad
  mito.upper.limit <- cur.med + diff.val

  get_quantile_perc_mt<-ecdf(seur@meta.data$percent.mt)
  cat("The percentage of cells will be removed based on mito.upper.limit is: ", 100*get_quantile_perc_mt(mito.upper.limit))
  #100*get_quantile_perc_mt(5)
  seur_fil <- subset(seur, subset = nFeature_RNA > gene.lower.limit & nFeature_RNA < gene.upper.limit &
                       nCount_RNA > umi.lower.limit & nCount_RNA < umi.upper.limit & percent.mt < mito.upper.limit)
  cell_filter_cutoff <- c(umi.upper.limit, umi.lower.limit, gene.upper.limit, gene.lower.limit, mito.upper.limit)
  return(list(cell_filter_cutoff, seur_fil))
}

