#!/usr/bin/env Rscript

# DM_cgmaptools.R for DMRichR
# Author: Shuye Pu
# Contributors: 

# Initialize --------------------------------------------------------------

cat("\n[DMRichR] Initializing \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")


suppressPackageStartupMessages(library(DMRichR))

# Global variables --------------------------------------------------------


cat("\n[DMRichR] Processing arguments from script \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

option_list <- list(
  optparse::make_option(c("-g", "--genome"), type = "character", default = NULL,
                        help = "Choose a genome (hg38, hg19, mm10, mm9, rheMac10, rheMac8, rn6, danRer11, galGal6, bosTau9, panTro6, dm6, susScr11, canFam3, TAIR10, or TAIR9) [required]"),
  optparse::make_option(c("-n", "--minSites"), type = "integer", default = 3,
                        help = "Choose the minimum number of Cytosines for a DMR [default = %default]"),
  optparse::make_option(c("-p", "--pval_cutoff"), type = "double", default = 0.05,
                        help = "Choose the cutoff value [from 0 to 1] for the pval of DML used for DMR detection [default = %default]"),
  optparse::make_option(c("-c", "--cores"), type = "integer", default = 20,
                        help = "Choose number of cores [default = %default]"),
  optparse::make_option(c("-r", "--ratio_cutoff"), type = "double", default = 2,
                        help = "Choose the cutoff value [from 0 to inf] for the ratio areaStat/nSites used for DMR detection [default = %default]"),
  optparse::make_option(c("-m", "--context"), type = "character", default = "CG",
                        help = "name of the methylation context [default = %default]"),
  optparse::make_option(c("-f", "--filepPatten"), type = "character", default = NULL,
                        help = "suffix of the report.txt.gz file [default = %default]"),
  optparse::make_option(c("-s", "--resPath"), type = "character", default = NULL,
                        help = "path to local resources [default = %default]"),
  optparse::make_option(c("-i", "--reportPath"), type = "character", default = NULL,
                        help = "path to the report.txt.gz file [default = %default]"),
  optparse::make_option(c("-o", "--wd"), type = "character", default = ".",
                        help = "path to the directory where DMRichR results will be stored [default = %default]"),
  optparse::make_option(c("-1", "--factor1"), type = "character", default = NULL,
                        help = "the first factor to be tested [default = %default]"),
  optparse::make_option(c("-2", "--factor2"), type = "character", default = NULL,
                        help = "the second factor to be tested in multifactor model [default = %default]")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# DM.R --------------------------------------------------------------------

DMRichR::DM_DSS.R(genome = opt$genome,
                 minSites =  opt$minSites,
                 pval_cutoff = opt$pval_cutoff,
                 ratio_cutoff = opt$raio_cutoff,
                 factor1 = opt$factor1,
                 factor2 = opt$factor2,
                 wd = opt$wd,
                 context = opt$context,
                 reportPath = opt$reportPath,
                 cores = opt$cores,
                 filePattern = opt$filePattern,
                 resPath = opt$resPath)
