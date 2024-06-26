% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DMRich.R
\name{DMRichCpG}
\alias{DMRichCpG}
\title{CpG annotation enrichment testing}
\usage{
DMRichCpG(
  sigRegions = sigRegions,
  regions = regions,
  genome = genome,
  resPath = NULL
)
}
\arguments{
\item{sigRegions}{\code{GRanges} object of DMRs}

\item{regions}{\code{GRanges} object of background regions}

\item{genome}{A character vector specifying the genome of interest
c("hg38", "hg19", "mm10", "mm9", "rheMac10", "rheMac8", "rn6", "danRer11",
 "galGal6", "bosTau9", "panTro6", "dm6", "susScr11", "canFam3")}

\item{resPath}{character specifying path to local resources in internet is not available}
}
\value{
A tibble with the enrichment results
}
\description{
Test DMRs for overlaps with CpG annotations (islands, shores, shelves, open sea)
}
\details{
DMRichCpG
}
