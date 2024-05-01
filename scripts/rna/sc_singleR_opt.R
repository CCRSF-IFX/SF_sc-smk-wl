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
setwd(opt$outdir)                        # %exclude_jupyterlab% 
saveRDS(opt, "opt.rds")                  # %exclude_jupyterlab%
