#!/usr/bin/env Rscript

# DM_cgmaptools.R for DMRichR
# Author: Shuye Pu
# Contributors: 

# Initialize --------------------------------------------------------------

cat("\n[DMRichR] Initializing \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

if(length(grep("genomecenter.ucdavis.edu", .libPaths())) > 0){
  .libPaths("/share/lasallelab/programs/DMRichR/R_4.1")
  AnnotationHub::setAnnotationHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_4.1")
  ExperimentHub::setExperimentHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_4.1")
}

if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
if(suppressPackageStartupMessages(!requireNamespace("DMRichR", quietly = TRUE))){
  Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = TRUE)
  remotes::install_github("shuye2009/DMRichR")
}
suppressPackageStartupMessages(library(DMRichR))

# Global variables --------------------------------------------------------

cat("\n[DMRichR] Processing arguments from script \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

option_list <- list(
  optparse::make_option(c("-g", "--genome"), type = "character", default = NULL,
                        help = "Choose a genome (hg38, hg19, mm10, mm9, rheMac10, rheMac8, rn6, danRer11, galGal6, bosTau9, panTro6, dm6, susScr11, canFam3, TAIR10, or TAIR9) [required]"),
  optparse::make_option(c("-n", "--minSites"), type = "integer", default = 5,
                        help = "Choose the minimum number of Cytosines for a DMR [default = %default]"),
  optparse::make_option(c("-o", "--cutoff"), type = "double", default = 0.05,
                        help = "Choose the cutoff value [from 0 to 1] for the single CpG coefficient utilized to discover testable background regions [default = %default]"),
  optparse::make_option(c("-c", "--cores"), type = "integer", default = 20,
                        help = "Choose number of cores [default = %default]"),
  optparse::make_option(c("-d", "--EnsDb"), type = "logical", default = FALSE,
                        help = "Logical to select Ensembl transcript annotation database [default = %default]"),
  optparse::make_option(c("-f", "--GOfuncR"), type = "logical", default = TRUE,
                        help = "Logical to run GOfuncR GO analysis [default = %default]"),
  optparse::make_option(c("--fileName"), type = "character", default = NULL,
                        help = "CGMaptools DMR file names [default = %default]"),
  optparse::make_option(c("--resPath"), type = "character", default = NULL,
                        help = "Path to local resources [default = %default]")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# DM.R --------------------------------------------------------------------

DMRichR::DM_cgmaptools.R(genome = opt$genome,
              minSites =  opt$minSites,
              cutoff = opt$cutoff,
              cores = opt$cores,
              GOfuncR = opt$GOfuncR,
              EnsDb = opt$EnsDb,
              fileName = opt$fileName,
              resPath = opt$resPath)
