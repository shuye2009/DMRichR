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
  optparse::make_option(c("-c", "--cores"), type = "integer", default = 10,
                        help = "Choose number of cores [default = %default]"),
  optparse::make_option(c("-r", "--ratio_cutoff"), type = "double", default = 2.0,
                        help = "Choose the cutoff value [from 0 to inf] for the ratio areaStat/nSites used for DMR detection [default = %default]"),
  optparse::make_option(c("-s", "--resPath"), type = "character", default = NULL,
                        help = "path to local resources [default = %default]"),
  optparse::make_option(c("-o", "--wd"), type = "character", default = ".",
                        help = "path to the directory where DMRichR results will be stored [default = %default]"),
  optparse::make_option(c("-1", "--factor1"), type = "character", default = "",
                        help = "the first factor to be tested [default = %default]"),
  optparse::make_option(c("-2", "--factor2"), type = "character", default = "",
                        help = "the second factor to be tested in multifactor model [default = %default]"),
  optparse::make_option(c("--ref1"), type = "character", default = "",
                        help = "reference condition for the first factor to be tested [default = %default]"),
  optparse::make_option(c("--ref2"), type = "character", default = "",
                        help = "reference condition for the second factor to be tested in multifactor model [default = %default]"),
  optparse::make_option(c("--condition1"), type = "character", default = "",
                        help = "reference condition for the two-group model [default = %default]"),
  optparse::make_option(c("--condition2"), type = "character", default = "",
                        help = "test condition for the two-group model [default = %default]"),
  optparse::make_option(c("--minDiff"), type = "double", default = 0.1,
                        help = "Choose the cutoff value [from 0 to 1] for the minimum methylation level difference between two groups [default = %default]"),
  optparse::make_option(c("--analysisType"), type = "character", default = "twoGroup",
                        help = "type of model of comparison, either 'twoGroup' or 'general' [default = %default]"),
  optparse::make_option(c("-d", "--override"), type = "logical", default = TRUE,
                        help = "whether to redefine DMR [default = %default]")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# DM.R --------------------------------------------------------------------

DMRichR::DSS.R(genome = opt$genome,
               analysisType = opt$analysisType,
               condition1 = opt$condition1,
               condition2 = opt$condition2,
               minDiff = opt$minDiff,
               minSites =  opt$minSites,
               pval_cutoff = opt$pval_cutoff,
               ratio_cutoff = opt$ratio_cutoff,
               factor1 = opt$factor1,
               factor2 = opt$factor2,
               ref1 = opt$ref1,
               ref2 = opt$ref2,
               wd = opt$wd,
               cores = opt$cores,
               resPath = opt$resPath,
               override = opt$override
             )
