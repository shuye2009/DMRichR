#' smoothPheatmap
#' @title DMR heatmap
#' @description Plot a heatmap of normalized individual smoothed methylation value z scores for selected regions (i.e. significant DMRs)
#' @param bs.filtered.bsseq Smoothed \code{bsseq} object
#' @param sigRegions \code{GRanges} object of regions to plot a heatmap for
#' @param testCovariate The factor tested for differences between groups
#' @param ... Additional arguments passed onto \link[pheatmap]{pheatmap}
#' @return Saves a pdf image of the heatmap in the DMR folder
#' @importFrom pheatmap pheatmap
#' @importFrom dplyr select_if
#' @importFrom RColorBrewer brewer.pal
#' @importFrom bsseq getMeth
#' @importClassesFrom bsseq BSseq 
#' @importMethodsFrom bsseq pData
#' @importFrom glue glue
#' @importFrom purrr pluck
#' @importFrom magrittr %>%
#' @references \url{https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/}
#' @export smoothPheatmap
#' 
smoothPheatmap <- function(bs.filtered.bsseq = bs.filtered.bsseq,
                           sigRegions = sigRegions,
                           testCovariate = testCovariate,
                           ...){
  cat("\n[DMRichR] DMR heatmap \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  bsseq::getMeth(BSseq = bs.filtered.bsseq,
                 regions = sigRegions,
                 type = "smooth",
                 what = "perRegion") %>% 
    na.omit() %>% 
    as.matrix() %>%
    pheatmap::pheatmap(.,
                       scale = "row",
                       annotation_col =  bs.filtered.bsseq %>% 
                         pData() %>%
                         as.data.frame() %>%
                         dplyr::select_if(~ nlevels(.) > 1),
                       color = RColorBrewer::brewer.pal(11,
                                                        name = "RdBu") %>%
                         rev(),
                       show_colnames = FALSE,
                       #angle_col = 45,
                       border_color = "grey",
                       main = glue::glue("Z-Scores of {length(sigRegions)} Differentially Methylated Regions"),
                       fontsize = 16,
                       filename = "./DMRs/heatmap.pdf",
                       width = 11,
                       height = 8.5,
                       annotation_colors = DMRichR::gg_color_hue(2) %>%
                         setNames(bs.filtered.bsseq %>%
                                    pData() %>%
                                    as.data.frame() %>%
                                    purrr::pluck(testCovariate) %>%
                                    unique() %>%
                                    sort() %>%
                                    rev()) %>% 
                                    list(testCovariate = .) %>% 
                         setNames(testCovariate),
                       ...
                       ) %>%
    return()
}
