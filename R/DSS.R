#' DSS.R
#' @title Run the pipeline with DSS as the DMR calling method
#' @description Performs the entire DMRichR analysis pipeline, 
#' which runs most functions in the package.
#' @param genome Character specifying the genome.
#' @param minDifff cutoff value [from 0 to inf] for the minimum difference between mean methylation 
#' levels between group1 and group2 during DMR detection for twoGroup analysis.
#' @param analysisType Character indicating type of analysis design, either 'general' or 'twoGroup'.
#' @param condition1 Character indicating the group1. 
#' @param condition2 Character indicating the group2.
#' @param pval_cutoff Numeric cutoff value [from 0 to 1] for the pval of DML used for DMR detection.
#' @param ratio_cutoff cutoff value [from 0 to inf] for the ratio areaStat/nSites used for DMR detection.
#' @param minSites Numeric for the minimum number of Cytosines for a DMR.
#' @param factor1 Character indicating factor of interest from the design matrix. 
#' @param factor2 Character indicating co-factor of interest from the design matrix.
#' @param ref1 Character indicating reference condition for the factor of interest from the design matrix. 
#' @param ref2 Character indicating reference condition for the co-factor of interest from the design matrix.
#' @param wd Character indicating the location where the analysis results to be stored.
#' @param cores Numeric specifying the number of cores to use. 10 is recommended. 
#' @param resPath character specifying path to local resources if internet is not available.
#' @param override logical indicating whether to redefine DMRs, default FALSE.
#' 
#' @importFrom dmrseq getAnnot dmrseq plotDMRs
#' @importFrom ggplot2 ggsave
#' @importFrom magrittr %>% %T>%
#' @importFrom purrr walk flatten set_names
#' @importFrom stringr str_remove str_trunc
#' @importFrom openxlsx read.xlsx write.xlsx
#' @importFrom BiocParallel MulticoreParam
#' @importFrom GenomicRanges GRangesList makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb dropSeqlevels seqlevels
#' @importFrom parallel mclapply
#' @importFrom glue glue
#' @importFrom dplyr mutate_if as_tibble pull case_when rename
#' @importFrom plyranges mutate filter
#' @importFrom forcats fct_rev
#' @importFrom Glimma glMDSPlot
#' @importFrom rGREAT submitGreatJob getEnrichmentTables plotRegionGeneAssociationGraphs
#' @importFrom enrichR listEnrichrDbs setEnrichrSite enrichr 
#' @importFrom utils write.table sessionInfo
#' @importFrom grDevices pdf dev.off
#' @importFrom parallel detectCores
#' @importFrom DSS DMLtest callDML callDMR
#' @importClassesFrom bsseq BSseq 
#' @importMethodsFrom bsseq pData seqnames sampleNames
#' @export DSS.R
#' 
DSS.R <- function(genome = c("hg38", "hg19", "mm10", "mm9", "rheMac10",
                            "rheMac8", "rn6", "danRer11", "galGal6",
                            "bosTau9", "panTro6", "dm6", "susScr11",
                            "canFam3", "TAIR10", "TAIR9"),
                 analysisType = "twoGroup",
                 condition1 = "",
                 condition2 = "",
                 minDiff = 0.1,
                 pval_cutoff = 0.05,
                 minSites = 3,
                 ratio_cutoff = 2,
                 factor1 = "",
                 factor2 = "",
                 ref1 = "",
                 ref2 = "",
                 cores = 10,
                 wd = ".",
                 resPath = NULL,
                 override = FALSE){
  
  if(!dir.exists(wd)) dir.create(wd)
  setwd(wd)
  # Check dmrseq version 
  if(Biobase::package.version("dmrseq") %>%
     stringr::str_remove("1.") %>%
     as.numeric() < 7.3){
    warning(paste("Your version of dmrseq is out of date and contains a bug.",
                  "This bug won't affect the DMRichR run but could affect your custom follow up analyses.",
                  "See the install section of the DMRichR README for the code to manually update.",
                  "Read more about the issue: https://github.com/kdkorthauer/dmrseq/issues/37"))
  }
  
  # Set options
  options(scipen = 999)
  options(readr.num_columns = 0)
  
  # Check for requirements
  stopifnot(genome %in% c("hg38", "hg19", "mm10", "mm9", "rheMac10",
                          "rheMac8", "rn6", "danRer11", "galGal6",
                          "bosTau9", "panTro6", "dm6", "susScr11",
                          "canFam3", "TAIR10", "TAIR9"))
  
  # Print
  print(glue::glue("genome = {genome}"))
  print(glue::glue("analysisType = {analysisType}"))
  print(glue::glue("minDiff = {minDiff}"))
  print(glue::glue("condition1 = {condition1}"))
  print(glue::glue("condition2 = {condition2}"))
  print(glue::glue("ratio_cutoff = {ratio_cutoff}"))
  print(glue::glue("pval_cutoff = {pval_cutoff}"))
  print(glue::glue("minSites = {minSites}"))
  print(glue::glue("factor1 = {factor1}"))
  print(glue::glue("factor2 = {factor2}"))
  print(glue::glue("ref1 = {ref1}"))
  print(glue::glue("ref2 = {ref2}"))
  print(glue::glue("cores = {cores}"))
  print(glue::glue("wd = {wd}"))
  print(glue::glue("resPath = {resPath}"))
  print(glue::glue("override = {override}"))
  
  # Experimental design ------------------------------------------------------
  if(analysisType == "general"){
    # factor1 and factor2 must be in the columns of design, must be releveled 
    design <- read.delim(file.path(wd, "sample_info.txt"), header = TRUE) |>
      dplyr::mutate(!!factor1:=factor(.data[[!!factor1]]), !!factor2:=factor(.data[[!!factor2]])) |>
      dplyr::mutate(!!factor1:=relevel(.data[[!!factor1]], ref1)) |>
      dplyr::mutate(!!factor2:=relevel(.data[[!!factor2]], ref2)) 
  }else if(analysisType == "twoGroup"){
    # factor1 and factor2 must be in the columns of design, must be releveled 
    design <- read.delim(file.path(wd, "sample_info.txt"), header = TRUE) |>
      dplyr::mutate(group=factor(group))
  }else{
    stop("analysis type ", analysisType, " is not supported!")
  }
  
  print(design)
  

  # Setup annotation databases ----------------------------------------------
  
  cat("\n[DMRichR] Selecting annotation databases \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  DMRichR::annotationDatabases(genome = genome,
                               EnsDb = FALSE)
  
  print(glue::glue("Saving Rdata..."))
  dir.create(file.path(wd, "RData"))
  settings_env <- ls(all = TRUE)
  save(list = settings_env, file = file.path(wd, "RData/settings.RData"))
  #load("RData/settings.RData")
  
  # Load and process samples ------------------------------------------------
  
  bs.filtered <- DMRichR::processReport(design, cores)
  
  print("pData")
  print(pData(bs.filtered))
  
  print(glue::glue("Building annotations for plotting..."))
  if(is(TxDb, "TxDb")){
    if(is.null(resPath)){
      annoTrack <- dmrseq::getAnnot(genome)
      saveRDS(annoTrack, file.path(wd, "RData/annoTrack.rds"))
    }else{
      annoTrack <- readRDS("~/resource/annoTrack.rds")
    }
    
  }else if(is(TxDb, "EnsDb")){
    annoTrack <- GenomicRanges::GRangesList(CpGs = DMRichR::getCpGs(genome, resPath),
                                            Exons = DMRichR::getExons(TxDb),
                                            compress = FALSE)
  }
  
  # Background --------------------------------------------------------------
  
  cat("\n[DMRichR] Getting bsseq background regions \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  dir.create("Extra")
  
  DMRichR::getBackground(bs.filtered,
                         minNumRegion = minSites,
                         maxGap = 1000) %>% 
    write.table(file = file.path(wd, "Extra/bsseq_background.csv"),
                sep = ",",
                quote = FALSE,
                row.names = FALSE)
  
  characterize_DMR <- function(DMR, dir){
    
    regions <- as(DMR$bgRegions, "GRanges")
    sigRegions <- as(DMR$sigRegions, "GRanges")
    
    print(glue::glue("Exporting DMR and background region information..."))
    dir.create("DMRs")
    output_DMR(DMR)
    
    if(sum(sigRegions$stat > 0) > 0 & sum(sigRegions$stat < 0) > 0){
      
      print(glue::glue("Summary: There are {tidySigRegions} DMRs \\
               ({tidyHyper}% hypermethylated, {tidyHypo}% hypomethylated) \\
               from {tidyRegions} background regions consisting of {tidyCpGs} CpGs \\
               assayed at 20x coverage", 
                       tidySigRegions = length(sigRegions),
                       tidyHyper = round(sum(sigRegions$stat > 0) / length(sigRegions), digits = 2)*100,
                       tidyHypo = round(sum(sigRegions$stat < 0) / length(sigRegions), digits = 2)*100,
                       tidyRegions = length(regions),
                       tidyCpGs = nrow(bs.filtered)))
    }
    
    
    # Annotate DMRs with gene symbols -----------------------------------------
    
    cat("\n[DMRichR] Annotating DMRs with gene symbols \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
    
    sigRegions %>%
      DMRichR::annotateRegions(TxDb = TxDb,
                               annoDb = annoDb,
                               resPath = resPath) %T>%
      openxlsx::write.xlsx(file = "DMRs/DMRs_annotated.xlsx")
    
    print(glue::glue("Annotating background regions with gene symbols..."))
    regions %>%
      DMRichR::annotateRegions(TxDb = TxDb,
                               annoDb = annoDb,
                               resPath = resPath) %>% 
      openxlsx::write.xlsx(file = "DMRs/background_annotated.xlsx")

    
    # HOMER -------------------------------------------------------------------
    
    sigRegions %>% 
      DMRichR::prepareHOMER(regions = regions, subfolder = dir)
    
    DMRichR::HOMER(genome = genome,
                   cores = cores,
                   subfolder = dir)
    
    # Smoothed global, chromosomal, and CGi methylation statistics ------------
    
    dir.create("Global")
    
    bs.filtered %>%
      DMRichR::globalStats(genome = genome,
                           testCovariate = NULL,
                           adjustCovariate = NULL,
                           matchCovariate = NULL,
                           resPath = resPath) %>%
      openxlsx::write.xlsx("Global/globalStats.xlsx")
    
    # Global plots ------------------------------------------------------------
    
    windows <- bs.filtered %>%
      DMRichR::windows(goi = goi)
    
    CpGs <- bs.filtered %>%
      DMRichR::CpGs()
    
    plots <- c("windows", "CpGs")
    
    if(genome %in% c("hg38", "hg19", "mm10", "mm9", "rheMac10", "rheMac8", "rn6", "danRer11", "galGal6",
                     "bosTau9", "panTro6", "dm6", "susScr11", "canFam3")){
      
      CGi <- bs.filtered %>% 
        DMRichR::CGi(genome = genome, resPath = resPath)
      
      plots <- c("windows", "CpGs", "CGi")
    }
    
    purrr::walk(plots,
                function(plotMatrix,
                         group =  bs.filtered %>%
                           pData() %>%
                           dplyr::as_tibble() %>%
                           dplyr::pull(!!factor1) %>%
                           forcats::fct_rev()){
                  
                  title <- dplyr::case_when(plotMatrix == "windows" ~ "20Kb Windows",
                                            plotMatrix == "CpGs" ~ "Single CpG",
                                            plotMatrix == "CGi" ~ "CpG Island")
                  
                  plotMatrix %>%
                    get() %>% 
                    DMRichR::PCA(testCovariate = factor1,
                                 bs.filtered = bs.filtered) %>%
                    ggplot2::ggsave(glue::glue("Global/{title} PCA.pdf"),
                                    plot = .,
                                    device = NULL,
                                    width = 11,
                                    height = 8.5)
                  
                  plotMatrix %>%
                    get() %>% 
                    DMRichR::densityPlot(group = group) %>% 
                    ggplot2::ggsave(glue::glue("Global/{title} Density Plot.pdf"),
                                    plot = .,
                                    device = NULL,
                                    width = 11,
                                    height = 4)
                  
                })
    
    
    
    # CpG and genic enrichment testing ----------------------------------------
    
    cat("\n[DMRichR] Performing DMRichments \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
    
    DMRich <- function(x){
      
      if(genome %in% c("hg38", "hg19", "mm10", "mm9", "rheMac10", "rheMac8", "rn6", "danRer11", "galGal6", "bosTau9", "panTro6", "dm6", "susScr11", "canFam3")){
        print(glue::glue("Running CpG annotation enrichments for {names(dmrList)[x]}"))
        dmrList[x] %>% 
          DMRichR::DMRichCpG(regions = regions,
                             genome = genome,
                             resPath = resPath) %T>%
          openxlsx::write.xlsx(file = glue::glue("DMRichments/{names(dmrList)[x]}_CpG_enrichments.xlsx")) %>% 
          DMRichR::DMRichPlot(type = "CpG") %>% 
          ggplot2::ggsave(glue::glue("DMRichments/{names(dmrList)[x]}_CpG_enrichments.pdf"),
                          plot = ., 
                          width = 4,
                          height = 3)
      }
      
      print(glue::glue("Running gene region annotation enrichments for {names(dmrList)[x]}"))
      dmrList[x] %>% 
        DMRichR::DMRichGenic(regions = regions,
                             TxDb = TxDb,
                             annoDb = annoDb,
                             resPath = resPath) %T>%
        openxlsx::write.xlsx(file = glue::glue("DMRichments/{names(dmrList)[x]}_genic_enrichments.xlsx")) %>% 
        DMRichR::DMRichPlot(type = "genic") %>% 
        ggplot2::ggsave(glue::glue("DMRichments/{names(dmrList)[x]}_genic_enrichments.pdf"),
                        plot = ., 
                        width = 4,
                        height = 4)
    }
    
    dmrList <- sigRegions %>% 
      DMRichR::dmrList()
    
    dir.create("DMRichments")
    
    purrr::walk(seq_along(dmrList),
                DMRich)
    
    purrr::walk(dplyr::case_when(genome %in% c("hg38", "hg19", "mm10", "mm9", "rn6") ~ c("CpG", "genic"),
                                 TRUE ~ "genic") %>%
                  unique(),
                function(type){
                  
                  print(glue::glue("Creating DMRich MultiPlots for {type} annotations"))
                  
                  DMRichR::DMparseR(direction =  c("All DMRs",
                                                   "Hypermethylated DMRs",
                                                   "Hypomethylated DMRs"),
                                    type = type) %>%
                    DMRichR::DMRichPlot(type = type,
                                        multi = TRUE) %>% 
                    ggplot2::ggsave(glue::glue("DMRichments/{type}_multi_plot.pdf"),
                                    plot = .,
                                    device = NULL,
                                    height = dplyr::case_when(type == "genic" ~ 5,
                                                              type == "CpG" ~ 3.5),
                                    width = 7)
                })
    
    # Overlap with human imprinted genes --------------------------------------
    
    cat("\n[DMRichR] Testing for imprinted gene enrichment \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
    
    dmrList <- sigRegions %>% 
      DMRichR::dmrList()
    
    sink("DMRs/human_imprinted_gene_overlaps.txt")
    
    purrr::walk(seq_along(dmrList),
                function(x){
                  print(glue::glue("Analyzing {names(dmrList)[x]}"))
                  
                  dmrList[x] %>%
                    DMRichR::imprintOverlap(regions = regions,
                                            TxDb = TxDb,
                                            annoDb = annoDb,
                                            resPath = resPath)
                })
    
    sink()
    
    # Manhattan plot ----------------------------------------------------------
    
    tryCatch({
      regions %>%
        DMRichR::annotateRegions(TxDb = TxDb,
                                 annoDb = annoDb,
                                 resPath = resPath) %>% 
        DMRichR::Manhattan()
    }, 
    error = function(error_condition) {
      print(glue::glue("Manhattan plot error"))
      #setwd("..")
    })
    
    # Gene Ontology, pathway analyses --------------------------------------------------
    
    cat("\n[DMRichR] Performing GREAT analyses \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
    dir.create("GREAT")
    
    dmr <- DMR$sigRegions
    
    gobp <- rGREAT::read_gmt(file.path(resPath, "c5.go.bp.v2024.1.Hs.symbols.gmt"),
                         from = "SYMBOL", to = "ENTREZ", orgdb = "org.Hs.eg.db")
    hallmark <- rGREAT::read_gmt(file.path(resPath, "h.all.v2024.1.Hs.symbols.gmt"),
                         from = "SYMBOL", to = "ENTREZ", orgdb = "org.Hs.eg.db")
    kegg <- rGREAT::read_gmt(file.path(resPath, "c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt"),
                         from = "SYMBOL", to = "ENTREZ", orgdb = "org.Hs.eg.db")
    reactome <- rGREAT::read_gmt(file.path(resPath, "c2.cp.reactome.v2024.1.Hs.symbols.gmt"),
                     from = "SYMBOL", to = "ENTREZ", orgdb = "org.Hs.eg.db")
    
    genesets <- list(gopb=gobp, hallmark=hallmark, kegg=kegg, reactome=reactome)
    for(status in c("hyper", "hypo")){
      gr <- as(dmr[dmr$status == status,], "GRanges")
      lapply(names(genesets), function(x){
        geneset <- genesets[[x]]
        tryCatch({
          GREAT_analysis(gr, 
                         geneset = geneset,
                         genesetName = x, 
                         padj_cutoff = 0.2, 
                         status = status, 
                         dname = "DMR", 
                         geneset_cutoff = 200, 
                         genome = genome)
        }, 
        error = function(error_condition) {
          print(glue::glue("GREAT for {x} error"))
        })
      })
    }
    
    
    # SUMMARY -------------------------------------------------------------------
    
    cat("\n[DMRichR] Summary \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
    
    print(glue::glue("Summary: There were {dmrLength} DMRs that covered {sigRegionPercent} of the genome. \\
                     The DMRs were identified from {backgroundLength } background regions that covered {regionPercent} of the genome.
                     {tidyHyper} of the DMRs were hypermethylated, and {tidyHypo} were hypomethylated. \\
                     The methylomes consisted of {tidyCpGs} CpGs.", 
                     dmrLength = sigRegions %>%
                       length() %>%
                       formatC(format = "d", big.mark = ","),
                     backgroundLength = regions %>%
                       length() %>%
                       formatC(format = "d", big.mark = ","),
                     tidyHyper = (sum(sigRegions$stat > 0) / length(sigRegions)) %>%
                       scales::percent(),
                     tidyHypo = (sum(sigRegions$stat < 0) / length(sigRegions)) %>%
                       scales::percent(),
                     tidyCpGs = nrow(bs.filtered) %>%
                       formatC(format = "d", big.mark = ","),
                     genomeSize = goi %>%
                       seqinfo() %>%
                       GenomeInfoDb::keepStandardChromosomes() %>%
                       as.data.frame() %>%
                       purrr::pluck("seqlengths") %>%
                       sum(),
                     dmrSize = sigRegions %>%
                       dplyr::as_tibble() %>%
                       purrr::pluck("width") %>%
                       sum(),
                     backgroundSize = regions  %>%
                       dplyr::as_tibble() %>%
                       purrr::pluck("width") %>%
                       sum(),
                     sigRegionPercent = (dmrSize/genomeSize) %>%
                       scales::percent(accuracy = 0.01),
                     regionPercent = (backgroundSize/genomeSize) %>%
                       scales::percent(accuracy = 0.01)
    ))
  }
  
  # DMRs --------------------------------------------------------------------
  
  cat("\n[DMRichR] Testing for DMRs with DSS \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  start_time <- Sys.time()
  
  if(analysisType == "general"){
    if(!file.exists("DMR_list.RDS") || override){
      DMR_lists <- DSS_multifactor(bs.filtered, design, factor1, factor2, 
                                    pval_cutoff, ratio_cutoff, minSites)
      DMR_lists <- purrr::compact(DMR_lists) # remove null elements
      saveRDS(DMR_lists, "DMR_list.RDS")
    }else{
      DMR_lists <- readRDS("DMR_list.RDS")
    }
    
    for(aname in names(DMR_lists)){
      
      dir <- file.path(wd, aname)
      dir.create(dir)
      setwd(dir)
      DMR <- DMR_lists[[aname]]
      print(glue::glue("Prossessing factor {aname}..."))
      characterize_DMR(DMR, dir)
    }
  }else if(analysisType=="twoGroup"){
    aname <- paste0(condition2, "_vs_", condition1)
    dir <- file.path(wd, aname)
    dir.create(dir)
    setwd(dir)
    
    DMR <- DSS_pairwise(bs.filtered, condition1, condition2, pval_cutoff, minDiff, minSites=3, cores=cores)
    print(glue::glue("Prossessing comparison {aname}..."))
    characterize_DMR(DMR, dir)
  }else{
    print(glue::glue("Analysis type {analysisType} is not supported!"))
  }
  
  print(glue::glue("DMR timing...in seconds"))
  end_time <- Sys.time()
  end_time - start_time
  
  # End ---------------------------------------------------------------------
  writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
  if(file.exists("Rplots.pdf")){file.remove("Rplots.pdf")}
   
  print(glue::glue("Done..."))
}
