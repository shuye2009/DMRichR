#' DM.R
#' @title Run the pipeline
#' @description Performs the entire DMRichR analysis pipeline, 
#' which runs most functions in the package.
#' @param genome Character specifying the genome.
#' @param coverage Numeric specifying the CpG coverage cutoff (1x recommended).
#' @param perGroup Numeric indicating percent of samples per 
#' a group to apply the CpG coverage cutoff to (from 0 to 1).
#' @param minCpGs Numeric for minimum number of CpGs for a DMR.
#' @param maxPerms Numeric indicating number of permutations for the DMR analysis.
#' @param maxBlockPerms Numeric indicating number of permutations for the block analysis.
#' @param cutoff Numeric indicating the cutoff value for the single CpG coefficient 
#' utilized to discover testable background regions. Values range from 0 to 1 and 
#' 0.05 (5 percent) is the default. If you get too many DMRs you should try 0.1 (10 percent).
#' @param testCovariate Character indicating factor of interest from the design matrix. 
#' Only case vs control is supported. 
#' @param adjustCovariate Character vector indicating discrete and continuous 
#' variables to adjust for based on the design matrix. Multiple variables can be provided. 
#' @param matchCovariate Character indicating the variable in the design matrix to block 
#' for when constructing permutations. Only a single variable can be provided and it cannot 
#' also be an adjustCovariate.
#' @param cores Numeric specifying the number of cores to use. 20 is recommended. 
#' @param GOfuncR Logical indicating whether to run a GOfuncR GO analysis.
#' @param sexCheck Logical indicating whether to confirm sex of each sample. 
#' This is highly recommended if your analysis has males and females 
#' and will also drop the sex chromosomes. You should also include the sex variable as an
#' adjustCovariate. 
#' @param EnsDb Logical indicating whether to to select Ensembl transcript annotation database.
#' This is recommended for non-model organisms. 
#' @param filePattern character indicating cytosine report file name pattern
#' @param resPath character specifying path to local resources if internet is not available
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
#' @importClassesFrom bsseq BSseq 
#' @importMethodsFrom bsseq pData seqnames sampleNames
#' @export DM.R
#' 
DSS.R <- function(genome = c("hg38", "hg19", "mm10", "mm9", "rheMac10",
                            "rheMac8", "rn6", "danRer11", "galGal6",
                            "bosTau9", "panTro6", "dm6", "susScr11",
                            "canFam3", "TAIR10", "TAIR9"),
                 pval_cutoff = 0.05,
                 minCpGs = 3,
                 ratio_cutoff = 2,
                 factor1 = "",
                 factor2 = "",
                 cores = 10,
                 wd = ".",
                 filePattern = ".CG_report.txt.gz",
                 context = "CG",
                 reportPath = ".",
                 resPath = NULL){
  
  wd <- file.path(wd, context)
  dir.create(wd)
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
  stopifnot(!is.null(testCovariate))
  stopifnot(coverage >= 1)
  
  # Check for more permutations than samples, 
  # factor1 and factor2 must be in the columns of design
  design <- read.delim("sample_info.txt", header = TRUE) 
  nSamples <- design %>%
    nrow()
  
  if(nSamples < maxPerms){
    print(glue::glue("Warning: You have requested {maxPerms} permutations for the DMR analysis, \\
                   which is more than the {nSamples} samples you have. \\
                   maxPerms will now be changed to {nSamples}."))
    maxPerms <- nSamples
  }
  
  if(nSamples < maxBlockPerms){
    print(glue::glue("Warning: You have requested {maxBlockPerms} permutations for the block analysis, \\
                   which is more than the {nSamples} samples you have. \\
                   maxBlockPerms will now be changed to {nSamples}."))
    maxBlockPerms <- nSamples
  }
  
  rm(nSamples)
  
  # Print
  print(glue::glue("genome = {genome}"))
  print(glue::glue("ratio_cutoff = {ratio_cutoff}"))
  print(glue::glue("pval_cutoff = {pval_cutoff}"))
  print(glue::glue("minCpGs = {minCpGs}"))
  print(glue::glue("factor1 = {factor1}"))
  print(glue::glue("factor2 = {factor2}"))
  print(glue::glue("cores = {cores}"))
  print(glue::glue("wd = {wd}"))
  print(glue::glue("reportPath = {reportPath}"))
  print(glue::glue("filePattern = {filePattern}"))
  print(glue::glue("resPath = {resPath}"))

  # Setup annotation databases ----------------------------------------------
  
  cat("\n[DMRichR] Selecting annotation databases \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  DMRichR::annotationDatabases(genome = genome,
                               EnsDb = EnsDb)
  
  print(glue::glue("Saving Rdata..."))
  dir.create("RData")
  settings_env <- ls(all = TRUE)
  save(list = settings_env, file = "RData/settings.RData")
  #load("RData/settings.RData")
  
  # Load and process samples ------------------------------------------------
  
  bs.filtered <- DMRichR::processReport(report_path, file_pattern, design, cores)
  

  print(glue::glue("Building annotations for plotting..."))
  if(is(TxDb, "TxDb")){
    if(is.null(resPath)){
      annoTrack <- dmrseq::getAnnot(genome)
      saveRDS(annoTrack, file.path(getwd(), "RData/annoTrack.rds"))
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
                         minNumRegion = minCpGs,
                         maxGap = 1000) %>% 
    write.table(file = "Extra/bsseq_background.csv",
                sep = ",",
                quote = FALSE,
                row.names = FALSE)
  
  # DMRs --------------------------------------------------------------------
  
  cat("\n[DMRichR] Testing for DMRs with DSS \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  start_time <- Sys.time()
  
  DMR_lists <- DSS_multi_factor(bs.filtered, desgin, factor1, factor2, pval_cutoff, ratio_cutoff)
  
  for(aname in names(DMR_lists)){
    dir.create(aname)
    DMR <- DMR_lists[[aname]]
    print(glue::glue("Prossessing factor {aname}..."))
    
    regions <- as(DMR$bgRegions, "GRanges")
    sigRegions <- as(DMR$sigRegions, "GRanges")
    
    print(glue::glue("Exporting DMR and background region information..."))
    output_DMR(DMR, file.path(wd, aname))
    
    if(sum(sigRegions$stat > 0) > 0 & sum(sigRegions$stat < 0) > 0){
      
      print(glue::glue("Summary: There are {tidySigRegions} DMRs \\
               ({tidyHyper}% hypermethylated, {tidyHypo}% hypomethylated) \\
               from {tidyRegions} background regions consisting of {tidyCpGs} CpGs \\
               assayed at {coverage}x coverage", 
                       tidySigRegions = length(sigRegions),
                       tidyHyper = round(sum(sigRegions$stat > 0) / length(sigRegions), digits = 2)*100,
                       tidyHypo = round(sum(sigRegions$stat < 0) / length(sigRegions), digits = 2)*100,
                       tidyRegions = length(regions),
                       tidyCpGs = nrow(bs.filtered)))
    }
    
    print(glue::glue("DMR timing..."))
    end_time <- Sys.time()
    end_time - start_time
    
    print(glue::glue("Annotating DMRs and plotting..."))
    
    pdf(file.path(wd, aname, "DMRs/DMRs.pdf"), height = 4, width = 8)
    tryCatch({
      DMRichR::plotDMRs2(bs.filtered,
                         regions = sigRegions,
                         testCovariate = NULL,
                         extend = (end(sigRegions) - start(sigRegions) + 1)*2,
                         addRegions = sigRegions,
                         annoTrack = annoTrack,
                         regionCol = "#FF00001A",
                         lwd = 2,
                         qval = FALSE,
                         stat = FALSE,
                         horizLegend = FALSE)
    },
    error = function(error_condition) {
      print(glue::glue("Warning: One (or more) of your DMRs can't be plotted, \\
                        try again later by manually loading R Data and subsetting sigRegions"))
    })
    dev.off()
    
    # Annotate DMRs with gene symbols -----------------------------------------
    
    cat("\n[DMRichR] Annotating DMRs with gene symbols \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
    
    sigRegions %>%
      DMRichR::annotateRegions(TxDb = TxDb,
                               annoDb = annoDb,
                               resPath = resPath) %T>%
      DMRichR::DMReport(regions = regions,
                        bs.filtered = bs.filtered,
                        coverage = 1,
                        name = "DMReport") %>% 
      openxlsx::write.xlsx(file = file.path(wd, aname, "DMRs/DMRs_annotated.xlsx"))
    
    print(glue::glue("Annotating background regions with gene symbols..."))
    regions %>%
      DMRichR::annotateRegions(TxDb = TxDb,
                               annoDb = annoDb,
                               resPath = resPath) %>% 
      openxlsx::write.xlsx(file = file.path(wd, aname, "DMRs/background_annotated.xlsx"))
    
    # Individual smoothed values ----------------------------------------------
    
    cat("\n[DMRichR] Smoothing individual methylation values \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
    start_time <- Sys.time()
    
  
    # ChromHMM and Reference Epigenomes ---------------------------------------
    
    if(length(grep("genomecenter.ucdavis.edu", .libPaths())) > 0 & genome == "hg38"){
      
      dir.create("LOLA")
      setwd("LOLA")
      
      dmrList <- sigRegions %>% 
        DMRichR::dmrList()
      
      LOLA <- function(x){
        
        dir.create(names(dmrList)[x])
        setwd(names(dmrList)[x])
        
        dmrList[x] %>%
          DMRichR::chromHMM(regions = regions,
                            cores = floor(cores/3)) %>% 
          DMRichR::chromHMM_heatmap()
        
        dmrList[x] %>%
          DMRichR::roadmap(regions = regions,
                           cores = floor(cores/3)) %>% 
          DMRichR::roadmap_heatmap()
        
        if(file.exists("Rplots.pdf")){file.remove("Rplots.pdf")}
      }
      
      parallel::mclapply(seq_along(dmrList),
                         LOLA,
                         mc.cores = 3,
                         mc.silent = TRUE)
      
      setwd("..")
    }
    
    # HOMER -------------------------------------------------------------------
    
    sigRegions %>% 
      DMRichR::prepareHOMER(regions = regions, subfolder = file.path(wd, aname))
    
    DMRichR::HOMER(genome = genome,
                   cores = cores,
                   subfolder = file.path(wd, aname))
    
    # Smoothed global, chromosomal, and CGi methylation statistics ------------
    
    dir.create("Global")
    
    bs.filtered %>%
      DMRichR::globalStats(genome = genome,
                           testCovariate = factor1,
                           adjustCovariate = factor2,
                           matchCovariate = NULL,
                           resPath = resPath) %>%
      openxlsx::write.xlsx(file.path(wd, aname, "Global/smoothed_globalStats.xlsx")) 
    
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
                           dplyr::pull(!!testCovariate) %>%
                           forcats::fct_rev()){
                  
                  title <- dplyr::case_when(plotMatrix == "windows" ~ "20Kb Windows",
                                            plotMatrix == "CpGs" ~ "Single CpG",
                                            plotMatrix == "CGi" ~ "CpG Island")
                  
                  plotMatrix %>%
                    get() %>% 
                    DMRichR::PCA(testCovariate = testCovariate,
                                 bs.filtered = bs.filtered) %>%
                    ggplot2::ggsave(glue::glue("{wd}/{aname}/Global/{title} PCA.pdf"),
                                    plot = .,
                                    device = NULL,
                                    width = 11,
                                    height = 8.5)
                  
                  plotMatrix %>%
                    get() %>% 
                    DMRichR::densityPlot(group = group) %>% 
                    ggplot2::ggsave(glue::glue("{wd}/{aname}/Global/{title} Density Plot.pdf"),
                                    plot = .,
                                    device = NULL,
                                    width = 11,
                                    height = 4)
                  
                  Glimma::glMDSPlot(plotMatrix %>%
                                      get(),
                                    groups = cbind(bsseq::sampleNames(bs.filtered),
                                                   pData(bs.filtered)) %>%
                                      dplyr::as_tibble() %>% 
                                      dplyr::select(-col) %>%
                                      dplyr::rename(Name = bsseq..sampleNames.bs.filtered.),
                                    path = file.path(wd, aname),
                                    folder = "interactiveMDS",
                                    html = glue::glue("{title} MDS plot"),
                                    launch = FALSE)
                })
    
    # Heatmap -----------------------------------------------------------------
    
    sigRegions %>%
      DMRichR::smoothPheatmap(bs.filtered = bs.filtered,
                              testCovariate = factor1)
    
    # CpG and genic enrichment testing ----------------------------------------
    
    cat("\n[DMRichR] Performing DMRichments \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
    
    DMRich <- function(x, dir){
      
      if(genome %in% c("hg38", "hg19", "mm10", "mm9", "rheMac10", "rheMac8", "rn6", "danRer11", "galGal6", "bosTau9", "panTro6", "dm6", "susScr11", "canFam3")){
        print(glue::glue("Running CpG annotation enrichments for {names(dmrList)[x]}"))
        dmrList[x] %>% 
          DMRichR::DMRichCpG(regions = regions,
                             genome = genome,
                             resPath = resPath) %T>%
          openxlsx::write.xlsx(file = glue::glue("{dir}/DMRichments/{names(dmrList)[x]}_CpG_enrichments.xlsx")) %>% 
          DMRichR::DMRichPlot(type = "CpG") %>% 
          ggplot2::ggsave(glue::glue("{dir}/DMRichments/{names(dmrList)[x]}_CpG_enrichments.pdf"),
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
        openxlsx::write.xlsx(file = glue::glue("{dir}/DMRichments/{names(dmrList)[x]}_genic_enrichments.xlsx")) %>% 
        DMRichR::DMRichPlot(type = "genic") %>% 
        ggplot2::ggsave(glue::glue("{dir}/DMRichments/{names(dmrList)[x]}_genic_enrichments.pdf"),
                        plot = ., 
                        width = 4,
                        height = 4)
    }
    
    dmrList <- sigRegions %>% 
      DMRichR::dmrList()
    
    dir.create(file.path(wd, aname, "DMRichments"))
    
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
                    ggplot2::ggsave(glue::glue("{dir}/DMRichments/{type}_multi_plot.pdf"),
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
    
    sink("{dir}/DMRs/human_imprinted_gene_overlaps.txt")
    
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
      setwd("..")
    })
    
    # Gene Ontology analyses --------------------------------------------------
    
    cat("\n[DMRichR] Performing gene ontology analyses \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
    
    dir.create("Ontologies")
    
    hyper <- sigRegions %>%
      plyranges::filter(stat > 0)
    hypo <- sigRegions %>%
      plyranges::filter(stat < 0)
    all_sigRegions <- list("hyper"=hyper, "hypo"=hypo)
    
    for(direction in c("hyper", "hypo")){
      asigRegions <- all_sigRegions[[direction]]
      
      if(genome %in% c("hg38", "hg19", "mm10", "mm9") & is.null(resPath)){
        
        print(glue::glue("Running GREAT"))
        GREATjob <- asigRegions %>%
          dplyr::as_tibble() %>%
          GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
          rGREAT::submitGreatJob(bg = regions,
                                 species = genome,
                                 rule = "oneClosest",
                                 request_interval = 1,
                                 version = "4.0.4")
        
        print(glue::glue("Saving and plotting GREAT results"))
        GREATjob %>%
          rGREAT::getEnrichmentTables(category = "GO") %T>% #%>%
          #purrr::map(~ dplyr::filter(., Hyper_Adjp_BH < 0.05)) %T>%
          openxlsx::write.xlsx(file = glue::glue("Ontologies/GREAT_results_{direction}.xlsx")) %>%
          DMRichR::slimGO(tool = "rGREAT",
                          annoDb = annoDb,
                          plots = FALSE) %T>%
          openxlsx::write.xlsx(file = glue::glue("Ontologies/GREAT_slimmed_results_{direction}.xlsx")) %>%
          DMRichR::GOplot() %>%
          ggplot2::ggsave(glue::glue("Ontologies/GREAT_plot_{direction}.pdf"),
                          plot = .,
                          device = NULL,
                          height = 8.5,
                          width = 10)
        
        # pdf(glue::glue("Ontologies/GREAT_gene_associations_graph.pdf"),
        #     height = 8.5,
        #     width = 11)
        # par(mfrow = c(1, 3))
        # res <- rGREAT::plotRegionGeneAssociationGraphs(GREATjob)
        # dev.off()
        # write.csv(as.data.frame(res),
        #           file = glue::glue("Ontologies/GREATannotations.csv"),
        #           row.names = FALSE)
        
      }
      
      if(GOfuncR == TRUE){
        print(glue::glue("Running GOfuncR"))
        asigRegions %>% 
          DMRichR::GOfuncR(regions = regions,
                           n_randsets = 1000,
                           upstream = 5000,
                           downstream = 1000,
                           annoDb = annoDb,
                           TxDb = TxDb) %T>%
          openxlsx::write.xlsx(glue::glue("Ontologies/GOfuncR_{direction}.xlsx")) %>% 
          DMRichR::slimGO(tool = "GOfuncR",
                          annoDb = annoDb,
                          plots = FALSE) %T>%
          openxlsx::write.xlsx(file = glue::glue("Ontologies/GOfuncR_slimmed_results_{direction}.xlsx")) %>% 
          DMRichR::GOplot() %>% 
          ggplot2::ggsave(glue::glue("Ontologies/GOfuncR_plot_{direction}.pdf"),
                          plot = .,
                          device = NULL,
                          height = 8.5,
                          width = 10)
      }
      
      
      if(genome != "TAIR10" & genome != "TAIR9" & is.null(resPath)){
        tryCatch({
          print(glue::glue("Running enrichR"))
          
          enrichR:::.onAttach() # Needed or else "EnrichR website not responding"
          #dbs <- enrichR::listEnrichrDbs()
          dbs <- c("GO_Biological_Process_2018",
                   "GO_Cellular_Component_2018",
                   "GO_Molecular_Function_2018",
                   "KEGG_2019_Human",
                   "Panther_2016",
                   "Reactome_2016",
                   "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")
          
          if(genome %in% c("mm10", "mm9", "rn6")){
            dbs %>%
              gsub(pattern = "Human", replacement = "Mouse")
          }else if(genome %in% c("danRer11", "dm6")){
            if(genome == "danRer11"){
              enrichR::setEnrichrSite("FishEnrichr")
            }else if(genome == "dm6"){
              enrichR::setEnrichrSite("FlyEnrichr")}
            dbs <- c("GO_Biological_Process_2018",
                     "GO_Cellular_Component_2018",
                     "GO_Molecular_Function_2018",
                     "KEGG_2019")
          }
          
          asigRegions %>%
            DMRichR::annotateRegions(TxDb = TxDb,
                                     annoDb = annoDb,
                                     resPath = resPath) %>%  
            dplyr::select(geneSymbol) %>%
            purrr::flatten() %>%
            enrichR::enrichr(dbs) %>% 
            purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis = "")) %T>% # %>% 
            #purrr::map(~ dplyr::filter(., Adjusted.P.value < 0.05)) %T>%
            openxlsx::write.xlsx(file = glue::glue("Ontologies/enrichr_{direction}.xlsx")) %>%
            DMRichR::slimGO(tool = "enrichR",
                            annoDb = annoDb,
                            plots = FALSE) %T>%
            openxlsx::write.xlsx(file = glue::glue("Ontologies/enrichr_slimmed_results_{direction}.xlsx")) %>% 
            DMRichR::GOplot() %>% 
            ggplot2::ggsave(glue::glue("Ontologies/enrichr_plot_{direction}.pdf"),
                            plot = .,
                            device = NULL,
                            height = 8.5,
                            width = 10)
          
        },
        error = function(error_condition) {
          print(glue::glue("Warning: enrichR did not finish. \\
                          The website may be down or there are internet connection issues."))
        })
      }
    }
    
    
    # Machine learning --------------------------------------------------------
    tryCatch({
      methylLearnOutput <- DMRichR::methylLearn(bs.filtered = bs.filtered,
                                                sigRegions = sigRegions,
                                                testCovariate = testCovariate,
                                                TxDb = TxDb,
                                                annoDb = annoDb,
                                                topPercent = 1,
                                                output = "all",
                                                saveHtmlReport = TRUE)
      
      if(!dir.exists("./Machine_learning")) {
        dir.create("./Machine_learning")
      } 
      
      if(length(methylLearnOutput) == 1) {
        openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput), 
                             file = "./Machine_learning/Machine_learning_output_one.xlsx") 
      } else {
        openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput$`Annotated common DMRs`,
                                  RF_Ranking_All_DMRs = methylLearnOutput$`RF ranking`,
                                  SVM_Ranking_All_DMRs = methylLearnOutput$`SVM ranking`),
                             file = "./Machine_learning/Machine_learning_output_all.xlsx") 
      }
      
      print(glue::glue("Saving RData..."))
      save(methylLearnOutput, file = "RData/machineLearning.RData")
      #load("RData/machineLearing.RData")
    },
    error = function(error_condition) {
      print(glue::glue("Warning: methylLearn did not finish. \\
                        You may have not had enough top DMRs across algrothims."))
    })
    
    # End ---------------------------------------------------------------------
    
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
    
    try(if(sum(blocks$pval < 0.05) > 0 & length(blocks) != 0){
      print(glue::glue("{length(sigBlocks)} significant blocks of differential methylation \\
             in {length(blocks)} background blocks"))
    }, silent = TRUE)
    
    writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
    if(file.exists("Rplots.pdf")){file.remove("Rplots.pdf")}
  } 
  print(glue::glue("Done..."))
}
