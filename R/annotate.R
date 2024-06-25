#' annotateRegions
#' @title Annotate DMRs and blocks
#' @description Annotate and tidy regions from \code{dmrseq::dmrseq()}
#' @param regions A \code{GRanges} object of DMRs, blocks, or background regions from \code{dmrseq::dmrseq()}
#' @param TxDb \code{TxDb} or \code{EnsDb} annotation package for genome of interest
#' @param annoDb Character specifying \code{OrgDb} annotation package for species of interest
#' @return A \code{tibble} of annotated regions
#' @rawNamespace import(ensembldb, except = c(select, filter))
#' @importFrom dplyr rename_with as_tibble case_when mutate select recode_factor distinct
#' @importFrom tidyselect any_of
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom ChIPseeker annotatePeak
#' @importFrom magrittr %>%
#' @importFrom glue glue glue_collapse
#' @importFrom GenomeInfoDb genome seqlevelsStyle
#' @importFrom plyranges join_overlap_left select mutate
#' @export annotateRegions
#' 
annotateRegions <- function(regions = sigRegions,
                            TxDb = TxDb,
                            annoDb = annoDb){
  
  genome <- TxDb %>%
    GenomeInfoDb::genome() %>%
    unique()
  
  print(glue::glue("Annotating {tidyRegions} regions from {genome} with gene symbols",
                   tidyRegions = length(regions)))
  
  if(is(TxDb, "EnsDb")){
    
    genome <- dplyr::case_when(GenomeInfoDb::genome(TxDb) == "GRCh38" ~ "hg38",
                               GenomeInfoDb::genome(TxDb) == "GRCm38" ~ "mm10",
                               GenomeInfoDb::genome(TxDb) == "Mmul_10" ~ "rheMac10",
                               GenomeInfoDb::genome(TxDb) == "Mmul_8.0.1" ~ "rheMac8",
                               GenomeInfoDb::genome(TxDb) == "Rnor_6.0" ~ "rn6",
                               GenomeInfoDb::genome(TxDb) == "GRCz11" ~ "danRer11",
                               GenomeInfoDb::genome(TxDb) == "GRCg6a" ~ "galGal6",
                               GenomeInfoDb::genome(TxDb) == "ARS-UCD1.2" ~ "bosTau9",
                               GenomeInfoDb::genome(TxDb) == "BDGP6.28" ~ "dm6",
                               GenomeInfoDb::genome(TxDb) == "Sscrofa11.1" ~ "susScr11",
                               GenomeInfoDb::genome(TxDb) == "CanFam3.1" ~ "canFam3") %>%
      unique()
                               }
  
  CpGs <- DMRichR::getCpGs(genome)
  
  regionsCpG <- regions %>% 
    plyranges::join_overlap_left(CpGs %>%
                                   plyranges::filter(type == "islands") %>% 
                                   plyranges::select(CpG.Island = type)) %>%
    unique() %>% 
    plyranges::join_overlap_left(CpGs %>%
                                   plyranges::filter(type == "shores") %>% 
                                   plyranges::select(CpG.Shore = type)) %>%
    unique() %>% 
    plyranges::join_overlap_left(CpGs %>%
                                   plyranges::filter(type == "shelves") %>% 
                                   plyranges::select(CpG.Shelf = type)) %>%
    unique() %>% 
    plyranges::join_overlap_left(CpGs %>%
                                   plyranges::filter(type == "inter") %>% 
                                   plyranges::select(Open.Sea = type)) %>%
    unique() %>% 
    plyranges::mutate(CpG.Island = dplyr::case_when(CpG.Island == "islands" ~ "Yes",
                                                    TRUE ~ "No"),
                      CpG.Shore = dplyr::case_when(CpG.Shore == "shores" ~ "Yes",
                                                   TRUE ~ "No"),
                      CpG.Shelf = dplyr::case_when(CpG.Shelf == "shelves" ~ "Yes",
                                                   TRUE ~ "No"),
                      Open.Sea = dplyr::case_when(Open.Sea == "inter" ~ "Yes",
                                                  TRUE ~ "No"))
  
  if(is(TxDb, "EnsDb")){
    GenomeInfoDb::seqlevelsStyle(regionsCpG) <- "Ensembl" # Work around for organism not supported
  }
  
  regionsCpG %>% 
    ChIPseeker::annotatePeak(TxDb = TxDb,
                             annoDb = annoDb,
                             overlap = "all",
                             verbose = FALSE) %>%
    dplyr::as_tibble() %>%
    dplyr::select(-tidyselect::any_of(c("strand",
                                        "index.start",
                                        "index.end",
                                        "index.width",
                                        "area",
                                        "geneChr",
                                        "geneStart",
                                        "geneEnd",
                                        "geneLength",
                                        "geneStrand",
                                        "transcriptId",
                                        "transcriptBiotype",
                                        "ENTREZID"))) %>%
    dplyr::mutate(annotation = gsub(" \\(.*","", annotation)) %>%
    dplyr::rename_with(
      ~ dplyr::case_when(
        . == "seqnames" ~ "chr",
        . == "L" ~ "CpGs",
        . == "beta" ~ "betaCoefficient",
        . == "stat" ~ "statistic",
        . == "pval" ~ "p.value",
        . == "qval" ~ "q.value",
        . == "SYMBOL" ~ "geneSymbol",
        . == "GENENAME" ~ "gene",
        TRUE ~ .)) %>% 
    return()
}

#' DMReport
#' @title Create an html report of DMRs or blocks
#' @description Create an html report of significant regions from \code{dmrseq}
#' @param sigRegions \code{GRanges} object of significant regions (DMRs or blocks) from \code{dmrseq} that 
#' were annotated by \code{DMRichR::annotateRegions}
#' @param regions \code{GRanges} object of background regions from \code{dmrseq}
#' @param bs.filtered Filtered \code{bsseq} object from \code{processBismark()}
#' @param coverage Numeric of coverage samples were filtered for
#' @param name Character for html report name
#' @return Saves an html report of DMRs with genic annotations
#' @importFrom gt gt tab_header fmt_number fmt_scientific fmt_percent as_raw_html
#' @importFrom dplyr select mutate 
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importClassesFrom bsseq BSseq 
#' @export DMReport
#' 
DMReport <- function(sigRegions = sigRegions,
                     regions = regions,
                     bs.filtered = bs.filtered,
                     coverage = coverage,
                     name = "DMReport"){
  cat("\n","Preparing HTML report...")
  
  stopifnot(class(sigRegions) == c("tbl_df", "tbl", "data.frame"))
  
  sigRegions %>%
    dplyr::select(chr,
                  start,
                  end,
                  width,
                  CpGs,
                  betaCoefficient,
                  statistic,
                  p.value,
                  q.value ,
                  difference,
                  CpG.Island,
                  CpG.Shore,
                  CpG.Shelf,
                  Open.Sea,
                  annotation,
                  distanceToTSS,
                  geneSymbol,
                  gene) %>%
    dplyr::mutate(difference = difference/100) %>% 
    gt::gt() %>%
    gt::tab_header(
      title = name,
      subtitle = glue::glue("There are {tidySigRegions} regions \\
             ({tidyHyper}% hypermethylated, {tidyHypo}% hypomethylated) \\
             from {tidyRegions} background regions consisting of {tidyCpGs} CpGs \\
             assayed at {coverage}x coverage.
             On average, the DMRs are {avgLength} bp long and contain {avgCpGs} CpGs.", 
                            tidySigRegions = nrow(sigRegions),
                            tidyHyper = round(sum(sigRegions$statistic > 0) / nrow(sigRegions),
                                              digits = 2)*100,
                            tidyHypo = round(sum(sigRegions$statistic < 0) / nrow(sigRegions),
                                             digits = 2)*100,
                            tidyRegions = length(regions),
                            tidyCpGs = nrow(bs.filtered),
                            avgLength = mean(sigRegions$width) %>% round(),
                            avgCpGs = mean(sigRegions$CpGs) %>% round()
                            )) %>% 
    gt::fmt_number(
      columns = gt::vars("width", "CpGs"),
      decimals = 0
      ) %>% 
    gt::fmt_scientific(
      columns = gt::vars("p.value", "q.value"),
      decimals = 2
      ) %>%
    gt::fmt_percent(
      columns = gt::vars("difference"),
      drop_trailing_zeros = TRUE
      ) %>% 
    gt::as_raw_html(inline_css = FALSE) %>%
    write(glue::glue("{name}.html"))
  cat("Done", "\n")
}

#' getExons
#' @title Obtain exons for plotting
#' @description Obtain exon annotations from a \code{ensDb} object and format for \code{plotDMRs()} 
#' @param TxDb A \code{ensDb} object
#' @return A \code{GRanges} object of annotated exons for every gene with a symbol in the genome.
#' @rawNamespace import(ensembldb, except = c(select, filter))
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importFrom BiocGenerics unlist
#' @importFrom plyranges mutate select filter
#' @importFrom GenomeInfoDb genome
#' @references Based on \code{annotatr::build_gene_annots()},
#'  see: \url{https://github.com/rcavalcante/annotatr/blob/master/R/build_annotations.R}
#' @export getExons
#' 
getExons <- function(TxDb = TxDb){
  
  stopifnot(is(TxDb, "EnsDb"))
  
  message('Building exons...')
  
  exons <- TxDb %>%
    ensembldb::cdsBy(by = "tx",
                     columns = c("tx_id", "gene_id", "symbol") #, # listColumns(TxDb)
                     # filter = GeneBiotypeFilter("protein_coding")
                     ) %>%
    BiocGenerics::unlist(use.names = FALSE) %>%
    plyranges::mutate(id = glue::glue("CDS:{seq_along(.)}"),
                      type = glue::glue("{unique(genome(TxDb))}_genes_cds")
                      ) %>%
    plyranges::select(id, tx_id, gene_id, symbol, type) # %>%
    # plyranges::filter(symbol != "")
  
  GenomeInfoDb::genome(exons) <- NA  
  ensembldb::seqlevelsStyle(exons) <- "UCSC"
  
  return(exons)
}

#' getCpGs
#' @title Obtain CpG island, CpG shore, CpG shelf, and open sea annotations
#' @description Obtain UCSC CpG islands and build CpG shore, CpG shelf, and open sea annotations.
#'  This function is based on \code{annotatr:::build_cpg_annots()}; however, 
#'  it obtains annotations for all genomes in the UCSC genome browser.
#' @param genome Character specifying the genome
#' @return A \code{GRanges} object of CpG island, CpG shore, CpG shelf, and open sea annotations
#' @importFrom GenomicRanges makeGRangesFromDataFrame mcols trim setdiff sort gaps
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importFrom plyranges stretch mutate select
#' @references Based on \code{annotatr:::build_cpg_annots()},
#'  see: \url{https://github.com/rcavalcante/annotatr/blob/master/R/build_annotations.R}
#' @export getCpGs
#' 
getCpGs <- function(genome = genome){
  
  message('Building CpG islands...')
  
  #islands <- readr::read_tsv(glue::glue("http://hgdownload.cse.ucsc.edu/goldenpath/{genome}/database/cpgIslandExt.txt.gz"),
  islands <- readr::read_tsv("extdata/cpgIslandExt.txt.gz",
                             col_names = c('chr','start','end'),
                             col_types = '-cii-------') %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
    GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
    plyranges::mutate(id = glue::glue("island:{seq_along(.)}"),
                      type = "islands")
  
  message('Building CpG shores...')
  
  shores <- islands %>% 
    plyranges::stretch(4000) %>% 
    GenomicRanges::trim() %>%
    GenomicRanges::setdiff(islands) %>%
    plyranges::mutate(id = glue::glue("shore:{seq_along(.)}"),
                      type = "shores")
  
  message('Building CpG shelves...')
  
  shelves <- shores %>% 
    plyranges::stretch(4000) %>% 
    GenomicRanges::trim() %>%
    GenomicRanges::setdiff(islands) %>%
    GenomicRanges::setdiff(shores) %>%
    plyranges::mutate(id = glue::glue("shelf:{seq_along(.)}"),
                      type = "shelves")
  
  message('Building inter-CpG-islands...')
  
  inter_cgi <- c(islands, shores, shelves) %>%
    GenomicRanges::sort() %>%
    GenomicRanges::gaps() %>%
    plyranges::mutate(id = glue::glue("inter:{seq_along(.)}"),
                      type = "inter")
  
  c(islands, shores, shelves, inter_cgi) %>%
    GenomicRanges::sort() %>%
    plyranges::mutate(tx_id = NA,
                      gene_id = NA,
                      symbol = NA) %>%
    plyranges::select(id, tx_id, gene_id, symbol, type) %>% 
    return()
}
