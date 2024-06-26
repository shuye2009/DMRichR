#' DM_cgmaptools.R
#' @title Run the pipeline for dmr derived from CGmaptools with single replicate
#' @description Performs the entire DMRichR analysis pipeline, 
#' using a limited number of tools.
#' @param genome Character specifying the genome.
#' @param minSites Numeric for minimum number of cytosine sites for a DMR.
#' @param cutoff Numeric cutoff for qvalues.
#' @param cores Numeric specifying the number of cores to use. 20 is recommended. 
#' @param GOfuncR Logical indicating whether to run a GOfuncR GO analysis.
#' @param fileName Character indicating dmr file name.
#' @param EnsDb Whether to use Ensemble database.
#' @param resPath character specifying path to local resources if internet is not available.
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
#' @export DM_cgmaptools.R
#' setwd("C:/PROJECTS/Shane/Harding_240124/local")
#' genome = "hg38"; minSites = 5; cutoff = 0.05; cores = 20; GOfuncR = TRUE; fileName = "WT_IR_vs_R172K_IR_dmr.CG.txt.gz"; 
#' EnsDb = FALSE; resPath = "C:/PROJECTS/Shane/Harding_240124/resource"

DM_cgmaptools.R <- function(genome = c("hg38", "hg19", "mm10", "mm9", "rheMac10",
                            "rheMac8", "rn6", "danRer11", "galGal6",
                            "bosTau9", "panTro6", "dm6", "susScr11",
                            "canFam3", "TAIR10", "TAIR9"),
                 minSites = 5,
                 cutoff = 0.05,
                 cores = 20,
                 GOfuncR = TRUE,
                 fileName = "*_dmr.CG.txt.gz",
                 EnsDb = FALSE,
                 resPath = NULL){
  

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
  print(glue::glue("minSites = {minSites}"))
  print(glue::glue("cutoff = {cutoff}"))
  print(glue::glue("cores = {cores}"))
  print(glue::glue("EnsDb = {EnsDb}"))
  print(glue::glue("GOfuncR = {GOfuncR}"))
  
  outfolder <- gsub(".txt.gz", "", fileName, fixed = TRUE)
  dir.create(outfolder)
  # Setup annotation databases ----------------------------------------------
  
  cat("\n[DMRichR] Selecting annotation databases \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  DMRichR::annotationDatabases(genome = genome,
                               EnsDb = EnsDb)
  
  print(glue::glue("Saving Rdata..."))
  dir.create(file.path(outfolder, "RData"))
  settings_env <- ls(all = TRUE)
  save(list = settings_env, file = file.path(outfolder, "RData", "settings.RData"))
  #load("RData/settings.RData")
  
  
  print(glue::glue("Building annotations for plotting..."))
  if(is(TxDb, "TxDb")){
    if(is.null(resPath)){
      annoTrack <- dmrseq::getAnnot(genome)
      saveRDS(annoTrack, file.path(outfolder, "RData","annoTrack.rds"))
    }else{
      annoTrack <- readRDS(file.path(resPath, "annoTrack.rds"))
    }
    
  }else if(is(TxDb, "EnsDb")){
    annoTrack <- GenomicRanges::GRangesList(CpGs = DMRichR::getCpGs(genome),
                                            Exons = DMRichR::getExons(TxDb),
                                            compress = FALSE)
  }

  
  # DMRs --------------------------------------------------------------------
  
  cat("\n[DMRichR] Processing DMRs from cgmaptools \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  start_time <- Sys.time()
  
  regions <- read.delim(fileName, header = FALSE) %>%
    `colnames<-`(c("chr", "start", "end", "stat", "pval", "beta", "pi", "nsites")) %>%
    dplyr::filter(!is.na(pval)) %>%
    dplyr::mutate(qval = p.adjust(pval)) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE, ignore.strand=TRUE)
  
  print(glue::glue("Selecting significant DMRs..."))
  
  regions <- regions %>% 
    plyranges::mutate(direction = dplyr::case_when(stat > 0 ~ "Hypermethylated",
                                                   stat < 0 ~ "Hypomethylated"),
                      difference = round((beta-pi)*100))
  
  if(sum(regions$qval < cutoff, na.rm=TRUE) < 100 & sum(regions$pval < cutoff, na.rm=TRUE) != 0){
    sigRegions <- regions %>%
      plyranges::filter(pval < cutoff)
  }else if(sum(regions$qval < cutoff, na.rm=TRUE) >= 100){
    sigRegions <- regions %>%
      plyranges::filter(qval < cutoff)
  }else if(sum(regions$pval < cutoff, na.rm=TRUE) == 0){
    stop(glue::glue("No significant DMRs detected in {length(regions)} background regions"))
  }

  
  print(glue::glue("Exporting DMR and background region information..."))
  
  dir.create(file.path(outfolder,"DMRs"))
  gr2bed(sigRegions, file.path(outfolder,"DMRs","DMRs.bed"))
  gr2bed(regions, file.path(outfolder,"DMRs", "backgroundRegions.bed"))
  
  if(sum(sigRegions$stat > 0) > 0 & sum(sigRegions$stat < 0) > 0){
    
    print(glue::glue("Summary: There are {tidySigRegions} DMRs \\
             ({tidyHyper}% hypermethylated, {tidyHypo}% hypomethylated) \\
             from {tidyRegions} background regions ", 
                     tidySigRegions = length(sigRegions),
                     tidyHyper = round(sum(sigRegions$stat > 0) / length(sigRegions), digits = 2)*100,
                     tidyHypo = round(sum(sigRegions$stat < 0) / length(sigRegions), digits = 2)*100,
                     tidyRegions = length(regions)))
  }
  
  print(glue::glue("DMR timing..."))
  end_time <- Sys.time()
  end_time - start_time
 

  
  # Annotate DMRs with gene symbols -----------------------------------------
  
  cat("\n[DMRichR] Annotating DMRs with gene symbols \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  sigRegions %>%
    DMRichR::annotateRegions(TxDb = TxDb,
                             annoDb = annoDb, resPath = resPath) %T>%
    openxlsx::write.xlsx(file = file.path(outfolder, "DMRs/DMRs_annotated.xlsx"))
  
  print(glue::glue("Annotating background regions with gene symbols..."))
  regions %>%
    DMRichR::annotateRegions(TxDb = TxDb,
                             annoDb = annoDb, resPath=resPath) %>% 
    openxlsx::write.xlsx(file = file.path(outfolder, "DMRs/background_annotated.xlsx"))
  

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
    DMRichR::prepareHOMER(regions = regions)
  
  DMRichR::HOMER(genome = genome,
                 cores = cores)
  system(paste0("mv HOMER ", outfolder))
  
  # CpG and genic enrichment testing ----------------------------------------
  
  cat("\n[DMRichR] Performing DMRichments \t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  DMRich <- function(x){
    
    if(genome %in% c("hg38", "hg19", "mm10", "mm9", "rheMac10", "rheMac8", "rn6", "danRer11", "galGal6", "bosTau9", "panTro6", "dm6", "susScr11", "canFam3")){
      print(glue::glue("Running CpG annotation enrichments for {names(dmrList)[x]}"))
      dmrList[x] %>% 
        DMRichR::DMRichCpG(regions = regions,
                           genome = genome, resPath = resPath) %T>%
        openxlsx::write.xlsx(file = glue::glue("{outfolder}/DMRichments/{names(dmrList)[x]}_CpG_enrichments.xlsx")) %>% 
        DMRichR::DMRichPlot(type = "CpG") %>% 
        ggplot2::ggsave(glue::glue("{outfolder}/DMRichments/{names(dmrList)[x]}_CpG_enrichments.pdf"),
                        plot = ., 
                        width = 4,
                        height = 3)
    }
    
    print(glue::glue("Running gene region annotation enrichments for {names(dmrList)[x]}"))
    dmrList[x] %>% 
      DMRichR::DMRichGenic(regions = regions,
                           TxDb = TxDb,
                           annoDb = annoDb) %T>%
      openxlsx::write.xlsx(file = glue::glue("{outfolder}/DMRichments/{names(dmrList)[x]}_genic_enrichments.xlsx")) %>% 
      DMRichR::DMRichPlot(type = "genic") %>% 
      ggplot2::ggsave(glue::glue("{outfolder}/DMRichments/{names(dmrList)[x]}_genic_enrichments.pdf"),
                      plot = ., 
                      width = 4,
                      height = 4)
  }
  
  dmrList <- sigRegions %>% 
    DMRichR::dmrList()
  
  dir.create(file.path(outfolder, "DMRichments"))
  
  purrr::walk(seq_along(dmrList),
              DMRich)
  
  purrr::walk(dplyr::case_when(genome %in% c("hg38", "hg19", "mm10", "mm9", "rn6") ~ c("CpG", "genic"),
                               TRUE ~ "genic") %>%
                unique(),
              function(type){
                
                print(glue::glue("Creating DMRich MultiPlots for {type} annotations"))
                
                DMRichR::DMparseR(subfolder = outfolder,
                                  direction =  c("All DMRs",
                                                 "Hypermethylated DMRs",
                                                 "Hypomethylated DMRs"),
                                  type = type) %>%
                  DMRichR::DMRichPlot(type = type,
                                      multi = TRUE) %>% 
                  ggplot2::ggsave(glue::glue("{outfolder}/DMRichments/{type}_multi_plot.pdf"),
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
  
  sink(file.path(outfolder, "DMRs/human_imprinted_gene_overlaps.txt"))
  
  purrr::walk(seq_along(dmrList),
              function(x){
                print(glue::glue("Analyzing {names(dmrList)[x]}"))
                
                dmrList[x] %>%
                  DMRichR::imprintOverlap(regions = regions,
                                          TxDb = TxDb,
                                          annoDb = annoDb)
              })
  
  sink()
  
  # Manhattan plot ----------------------------------------------------------
  wd <- getwd()
  tryCatch({
    regions %>%
      DMRichR::annotateRegions(TxDb = TxDb,
                               annoDb = annoDb, resPath = resPath) %>% 
      DMRichR::Manhattan(subfoler = outfolder)
  }, 
  error = function(error_condition) {
    print(glue::glue("Manhattan plot error"))
    setwd(wd)
  })
  
  # Gene Ontology analyses --------------------------------------------------
  
  cat("\n[DMRichR] Performing gene ontology analyses \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  dir.create(file.path(outfolder, "Ontologies"))
  
  if(genome %in% c("hg38", "hg19", "mm10", "mm9") & in.null(resPath)){
    
    print(glue::glue("Running GREAT"))
    GREATjob <- sigRegions %>%
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
      openxlsx::write.xlsx(file = glue::glue("Ontologies/GREAT_results.xlsx")) %>%
      DMRichR::slimGO(tool = "rGREAT",
                      annoDb = annoDb,
                      plots = FALSE) %T>%
      openxlsx::write.xlsx(file = glue::glue("Ontologies/GREAT_slimmed_results.xlsx")) %>%
      DMRichR::GOplot() %>%
      ggplot2::ggsave(glue::glue("Ontologies/GREAT_plot.pdf"),
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
    sigRegions %>% 
      DMRichR::GOfuncR(regions = regions,
                       n_randsets = 1000,
                       upstream = 5000,
                       downstream = 1000,
                       annoDb = annoDb,
                       TxDb = TxDb) %T>%
      openxlsx::write.xlsx(glue::glue("{outfolder}/Ontologies/GOfuncR.xlsx")) %>% 
      DMRichR::slimGO(tool = "GOfuncR",
                      annoDb = annoDb,
                      plots = FALSE) %T>%
      openxlsx::write.xlsx(file = glue::glue("{outfolder}/Ontologies/GOfuncR_slimmed_results.xlsx")) %>% 
      DMRichR::GOplot() %>% 
      ggplot2::ggsave(glue::glue("{outfolder}/Ontologies/GOfuncR_plot.pdf"),
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
      
      sigRegions %>%
        DMRichR::annotateRegions(TxDb = TxDb,
                                 annoDb = annoDb) %>%  
        dplyr::select(geneSymbol) %>%
        purrr::flatten() %>%
        enrichR::enrichr(dbs) %>% 
        purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis = "")) %T>% # %>% 
        #purrr::map(~ dplyr::filter(., Adjusted.P.value < 0.05)) %T>%
        openxlsx::write.xlsx(file = glue::glue("Ontologies/enrichr.xlsx")) %>%
        DMRichR::slimGO(tool = "enrichR",
                        annoDb = annoDb,
                        plots = FALSE) %T>%
        openxlsx::write.xlsx(file = glue::glue("Ontologies/enrichr_slimmed_results.xlsx")) %>% 
        DMRichR::GOplot() %>% 
        ggplot2::ggsave(glue::glue("Ontologies/enrichr_plot.pdf"),
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
  
 
  # End ---------------------------------------------------------------------
  
  cat("\n[DMRichR] Summary \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  
  print(glue::glue("Summary: There were {dmrLength} DMRs that covered {sigRegionPercent} of the genome. \\
                   The DMRs were identified from {backgroundLength } background regions that covered {regionPercent} of the genome.
                   {tidyHyper} of the DMRs were hypermethylated, and {tidyHypo} were hypomethylated.", 
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
 
  writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
  if(file.exists("Rplots.pdf")){file.remove("Rplots.pdf")}
  
  print(glue::glue("Done..."))
}
