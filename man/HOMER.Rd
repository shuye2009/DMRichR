% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HOMER.R
\name{HOMER}
\alias{HOMER}
\title{Transcription factor motif analysis}
\usage{
HOMER(genome = genome, cores = cores, subfolder = ".")
}
\arguments{
\item{genome}{Character specifying the genome}

\item{cores}{Integer specifying the number of cores to use}

\item{subfolder}{Character specifying saving directory to allow multiple instance of HOMER}
}
\value{
A folder with HOMER results
}
\description{
Run a HOMER known transcription motif analysis for all DMRs,
hypermethylated DMRs, and hypomethylated DMRs. The function requires HOMER's 
\code{findMotifsGenome.pl} script to be in the path (i.e. \code{module load homer})
and for the genome of interest (see \code{perl /path-to-homer/configureHomer.pl -list}) 
to be installed through \code{perl /path-to-homer/configureHomer.pl -install human}.
}
\details{
HOMER
}
\references{
\url{http://homer.ucsd.edu/homer/introduction/configure.html}

\url{https://www.biostars.org/p/443759/}
}
