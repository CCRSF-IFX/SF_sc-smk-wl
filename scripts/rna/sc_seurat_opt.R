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
setwd(opt$outdir)                        # %exclude_jupyterlab% 
saveRDS(opt, "opt.rds")                  # %exclude_jupyterlab%