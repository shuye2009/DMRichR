#' windows
#' @title Extract methylation values from tiled genomic windows
#' @description Obtain windows of individual smoothed methylation values
#' @param bs.filtered.bsseq Smoothed \code{bsseq} object with a testCovariate in \code{pData}
#' @param size The number of bases in the window (default is 2e4, which is 20 Kb)
#' @param goi A \code{BSgenome} object of the genome of interest (i.e. "BSgenome.Hsapiens.UCSC.hg38")
#' @return A matrix of smoothed individual methylation values
#' @importFrom magrittr %>%
#' @importFrom GenomeInfoDb seqlengths keepStandardChromosomes
#' @importFrom GenomicRanges tileGenome
#' @importFrom bsseq getMeth
#' @importClassesFrom bsseq BSseq 
#' @importMethodsFrom bsseq pData
#' @export windows
#' 
windows <- function(bs.filtered.bsseq = bs.filtered.bsseq,
                    size = 2e4,
                    goi = goi){
  print(glue::glue("Obtaining {size/1000} Kb window individual smoothed methylation values from the {BSgenome::commonName(goi)} genome"))
  goi %>%
    GenomeInfoDb::seqlengths() %>%
    GenomicRanges::tileGenome(tilewidth = size,
                              cut.last.tile.in.chrom = TRUE) %>%
    GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
    bsseq::getMeth(BSseq = bs.filtered.bsseq,
                   regions = .,
                   type = "raw",
                   what = "perRegion") %>% 
    na.omit() %>%
    return()
}

#' CGi
#' @title Extract methylation values from CpG islands
#' @description Obtain individual smoothed methylation values for CpG islands
#' @param bs.filtered.bsseq Smoothed \code{bsseq} object with a testCovariate in \code{pData}
#' @param genome A character vector of the genome of interest (i.e. "hg38")
#' @return A matrix of smoothed individual methylation values
#' @importFrom plyranges filter
#' @importFrom magrittr %>%
#' @importFrom bsseq getMeth
#' @importClassesFrom bsseq BSseq 
#' @importMethodsFrom bsseq pData
#' @export CGi
#' 
CGi <- function(bs.filtered.bsseq = bs.filtered.bsseq,
                genome = genome, resPath = resPath){
  
  print(glue::glue("Obtaining individual smoothed methylation values of CpG islands from {genome}"))
  
  genome %>%
    DMRichR::getCpGs(resPath = resPath) %>% 
    plyranges::filter(type == "islands") %>% 
    bsseq::getMeth(BSseq = bs.filtered.bsseq,
                   regions = .,
                   type = "raw",
                   what = "perRegion") %>% 
    na.omit() %>%
    return()
}

#' CpGs
#' @title Extract single CpG methylation values
#' @description Extracts single CpG individual smoothed methylation values
#' @param bs.filtered.bsseq Smoothed \code{bsseq} object with a testCovariate in \code{pData}
#' @return A matrix of smoothed individual methylation values
#' @importFrom bsseq getMeth
#' @importClassesFrom bsseq BSseq 
#' @importMethodsFrom bsseq pData
#' @export CpGs
#' 
CpGs <- function(bs.filtered.bsseq = bs.filtered.bsseq){
  
  print(glue::glue("Obtaining smoothed methylation values for all covered CpGs"))
  
  bs.filtered.bsseq %>% 
    bsseq::getMeth(BSseq = .,
                   type = "raw",
                   what = "perBase") %>%
    na.omit() %>%
    return()
}

#' PCA
#' @title PCA plot of extracted methylation values
#' @description Performs and plots a PCA from individual smoothed methylation values
#' @param matrix A matrix of smoothed individual methylation values
#' @param testCovariate Factor of interest
#' @param bs.filtered.bsseq Smoothed \code{bsseq} object with a testCovariate in \code{pData}
#' @return A \code{ggplot} object that can be viewed by calling it,
#'  saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax
#' @importFrom PCAtools pca biplot
#' @importFrom glue glue
#' @importFrom Glimma glMDSPlot
#' @importFrom purrr pluck
#' @importClassesFrom bsseq BSseq 
#' @importMethodsFrom bsseq pData
#' @references \url{https://stackoverflow.com/questions/40315227/how-to-solve-prcomp-default-cannot-rescale-a-constant-zero-column-to-unit-var/40317343}
#' @export PCA
#' 
PCA <- function(matrix = matrix,
                testCovariate = testCovariate,
                bs.filtered.bsseq = bs.filtered.bsseq){
  
  print(glue::glue("PCA of {length(matrix)} sites"))
  
  matrix %>%
    PCAtools::pca(scale = TRUE,
                  metadata = pData(bs.filtered.bsseq),
                  removeVar = 0.1) %>% 
    PCAtools::biplot(colby = testCovariate,
                     colkey = DMRichR::gg_color_hue(2) %>%
                       setNames(bs.filtered.bsseq %>%
                                  pData() %>%
                                  as.data.frame() %>%
                                  purrr::pluck(testCovariate) %>%
                                  unique() %>%
                                  sort() %>%
                                  rev()),
                     pointSize = 6,
                     labSize = 4,
                     legendPosition = 'top',
                     legendLabSize = 16,
                     legendIconSize = 8.0)
}

#' densityPlot
#' @title Density plot of extracted methylation values
#' @description Creates a density plot of the mean of individual smoothed methylation values
#' @param matrix A matrix of smoothed individual methylation values
#' @param group Ordered factor vector of sample groupings
#' @return A \code{ggplot} object that can be viewed by calling it,
#'  saved with \code{ggplot2::ggsave()}, or further modified by adding \code{ggplot2} syntax.
#' @importFrom magrittr %>% set_colnames
#' @importFrom dplyr as_tibble select transmute contains mutate
#' @importFrom forcats fct_rev
#' @import ggplot2
#' @importFrom tidyr gather
#' @export densityPlot
#' 
densityPlot <- function(matrix = matrix,
                        group = NA){
  
  print(glue::glue("Density plot of {length(matrix)} sites"))
  
  matrix  %>%
    dplyr::as_tibble() %>% 
    magrittr::set_colnames(paste(group, seq_along(1:length(group)))) %>%
    dplyr::transmute(Group1 = dplyr::select(., dplyr::contains(levels(group)[1])) %>% rowMeans()*100,
                     Group2 = dplyr::select(., dplyr::contains(levels(group)[2])) %>% rowMeans()*100) %>%
    magrittr::set_colnames(c(levels(group)[1], levels(group)[2])) %>% 
    tidyr::gather(key = "variable",
                  value = "value") %>%
    #dplyr::mutate(variable = factor(.$variable)) %>% 
    dplyr::mutate(variable = factor(.$variable, levels = levels(group))) %>% 
    ggplot(aes(value, color = variable)) +
    geom_density(size = 1.2) +
    labs(x = "Percent Methylation",
         y = "Density",
         color = "Group") +
    theme_classic() +
    scale_x_continuous(expand = c(0.05,0.05),
                       breaks = c(0,25,50,75,100)) +
    scale_y_continuous(expand = c(0.00,0.001)) +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          strip.text = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.position = "bottom",
          legend.title = element_text(size = 20)) %>%
    
    return()
}
