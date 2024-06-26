% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/globalPlots.R
\name{PCA}
\alias{PCA}
\title{PCA plot of extracted methylation values}
\usage{
PCA(
  matrix = matrix,
  testCovariate = testCovariate,
  bs.filtered.bsseq = bs.filtered.bsseq
)
}
\arguments{
\item{matrix}{A matrix of smoothed individual methylation values}

\item{testCovariate}{Factor of interest}

\item{bs.filtered.bsseq}{Smoothed \code{bsseq} object with a testCovariate in \code{pData}}
}
\value{
A \code{ggplot} object that can be viewed by calling it,
 saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax
}
\description{
Performs and plots a PCA from individual smoothed methylation values
}
\details{
PCA
}
\references{
\url{https://stackoverflow.com/questions/40315227/how-to-solve-prcomp-default-cannot-rescale-a-constant-zero-column-to-unit-var/40317343}
}
