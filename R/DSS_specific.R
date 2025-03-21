#' processReport
#' @title Process methylation report files
#' @description Build bsseq object from report files and perform smoothing
#' @param design sample information dataframe.
#' @param cores Numeric specifying the number of cores to use. 10 is recommended.
#' 
#' @return a bsseq object
#' @export processReport

processReport <- function(design, cores){
  files <- design$path
  if(Sys.info()['sysname'] == "Windows"){
    bpparams <- BiocParallel::SnowParam(workers = parallel::detectCores()-2, progressbar = TRUE)
  }else{
    bpparams <-  BiocParallel::MulticoreParam(workers = cores, progressbar = TRUE)
  }
  
  print(glue::glue("Reading cytosine reports..."))
  if(!file.exists("bs.RDS")){
    bs <- bsseq::read.bismark(files,
                              colData = design,
                              rmZeroCov = FALSE,
                              strandCollapse = TRUE,
                              verbose = TRUE,
                              BPPARAM = bpparams,
                              nThread = 1) # 1L # nThread
    bs <- GenomeInfoDb::keepStandardChromosomes(bs, pruning.mode = "coarse")
    GenomeInfoDb::seqlevelsStyle(bs) <- "UCSC"
    
    saveRDS(bs, "bs.RDS")
  }else{
    bs <- readRDS("bs.RDS")
  }
  
  return(bs)
}


#' DSS_multifactor
#' @title DML and DMR calling with DSS method
#' @description Call DML using multifactor design, then get DMRS for each factor 
#' and interaction effect if two factors are provided in the design
#' @param bss a smoothed bsseq object.
#' @param design sample information dataframe.
#' @param pval_cutoff Numeric cutoff value [from 0 to 1] for the pval of DML used for DMR detection.
#' @param ratio_cutoff cutoff value [from 0 to inf] for the ratio areaStat/nSites used for DMR detection.
#' @param minSites Numeric for the minimum number of Cytosines for a DMR.
#' @param factor1 Character indicating factor of interest from the design matrix. 
#' @param factor2 Character indicating co-factor of interest from the design matrix.
#' 
#' @return a list of DMRs for all factors and interaction. See \code{findDMR} for details
#' @export DSS_multifactor

DSS_multifactor <- function(bss, design, factor1, factor2, pval_cutoff, ratio_cutoff, minSites=3){
  if(!file.exists("DMLfit.RDS")){
    if(is.null(factor2)){
      fomu <- as.formula(paste("~", factor1))
    }else{
      fomu <- as.formula(paste("~", factor1, "+", factor2, "+", factor1, ":", factor2))
    }
    DMLfit= DSS::DMLfit.multiFactor(bss, 
                                    design = design, 
                                    formula = fomu, 
                                    smoothing = TRUE,
                                    smoothing.span = 500)
    saveRDS(DMLfit, "DMLfit.RDS")
  }else{
    DMLfit <- readRDS("DMLfit.RDS")
  }
  
  message("design matrix:")
  print(DMLfit$X)                
  
  message("[DSS_multi_factor] finding DML for factor ", factor1)
  DMLfactor1 <- DSS::DMLtest.multiFactor(DMLfit, term = factor1)
  if(!is.null(factor2)){
    message("[DSS_multi_factor] finding DML for factor ", factor2)
    DMLfactor2 <- DSS::DMLtest.multiFactor(DMLfit, term = factor2)
    message("[DSS_multi_factor] finding DML for interaction of ", factor1, " and ", factor2)
    DMLInter <- DSS::DMLtest.multiFactor(DMLfit, term = factor(paste0(factor1, ":", factor2)))
  }

  rm(DMLfit)
  
  message("[DSS_multi_factor] finding DMR for factor ", factor1)
  
  dmrsfactor1 <- findDMR(DMLfactor1, 
                        pval_cutoff = pval_cutoff, 
                        ratio_cutoff = ratio_cutoff,
                        minSites = minSites)
  if(!is.null(factor2)){
    message("[DSS_multi_factor] finding DMR for factor ", factor2)
    dmrsfactor2 <- findDMR(DMLfactor2, 
                          pval_cutoff = pval_cutoff, 
                          ratio_cutoff = ratio_cutoff,
                          minSites = minSites)
    message("[DSS_multi_factor] finding DMR for interaction of ", factor1, " and ", factor2)
    dmrsInter <- findDMR(DMLInter, 
                         pval_cutoff = pval_cutoff, 
                         ratio_cutoff = ratio_cutoff,
                         minSites = minSites)
  }else{
    dmrsfactor2 <- dmrsInter <- NULL
  }
  
  return(list("factor1"= dmrsfactor1, "factor2"= dmrsfactor2, "interaction" = dmrsInter))
}

#' findDMR
#' @title DMR calling from DML
#' @description Call DMR using multifactor design
#' @param DML a smoothed bsseq object.
#' @param pval_cutoff Numeric cutoff value [from 0 to 1] for the pval of DML used for DMR detection.
#' @param ratio_cutoff cutoff value [from 0 to inf] for the ratio areaStat/nSites used for DMR detection.
#' @param minSites Numeric for the minimum number of Cytosines for a DMR.
#' 
#' @return a list of two dataframes, the first is for signigicant DMRs and the second is for background DMRs
#' @export findDMR
#' 
findDMR <- function(DML, pval_cutoff=0.05, ratio_cutoff=2, minSites=3){
  
  DML <- DML |>
    dplyr::filter(!is.na(pvals))
  
  message("[findDMR] sample of significant DML at ", pval_cutoff, "...")
  print(head(DML[DML$pvals<pval_cutoff,]))
  
  dmrs <- DSS::callDMR(DML, 
                  delta=0, 
                  p.threshold=pval_cutoff,
                  minlen=50, 
                  minCG=minSites, 
                  dis.merge=100, 
                  pct.sig=0.5) |>
    dplyr::filter(!is.na(chr)) |>
    dplyr::mutate(ratio = areaStat/nCG) |>
    dplyr::filter(abs(ratio) > ratio_cutoff) |>
    dplyr::mutate(stat = areaStat, .keep = "unused") |>
    dplyr::mutate(status = case_when(stat > 0 ~ "hyper",
                                     stat < 0 ~ "hypo",
                                     .default = "none"))
  
  message("[findDMR] sample of significant DMR at pval_cutoff ", pval_cutoff, 
          " before filtering for ratio_cutoff ", ratio_cutoff)
  print(head(dmrs))
  
  if(nrow(dmrs) > 0){
    
    message("select background DMR")
    dmrs_background <- DSS::callDMR(DML, 
                               delta=0, 
                               p.threshold=1,
                               minlen=50, 
                               minCG=minSites, 
                               dis.merge=0, # do not merge regions
                               pct.sig=0.5) |>
      dplyr::filter(!is.na(chr)) |>
      dplyr::mutate(ratio = areaStat/nCG) |>
      dplyr::mutate(stat = areaStat, .keep = "unused") |>
      dplyr::mutate(status = case_when(stat > 0 ~ "hyper",
                                       stat < 0 ~ "hypo",
                                       .default = "none")) 
    
    return(list(sigRegions=dmrs, bgRegions=dmrs_background))
  }else{
    message("!!! No significant DMR found at p.threshold ", pval_cutoff)
    return(NULL)
  }
}

#' GREAT_analysis
#' @title GREAT analysis
#' @description Perform GREAT analysis for DMR
#' @param gr a GRanges object.
#' @param geneset a list from rGREAT::read_gmt.
#' @param padj_cutoff Numeric cutoff value [from 0 to 1] for the adjusted p-value.
#' @param geneset_cutoff cutoff value for the maximal size of the geneset.
#' @param status indicates the direction of methylation changes, either hyper or hypo.
#' @param dname indicates nature of DMR, either sigDMR or background
#' @param dir path for saving GREAT analysis results.
#' @param genome Character specifying the genome.
#' 
#' @export GREAT_analysis

GREAT_analysis <- function(gr, geneset, genesetName="GO:BP", padj_cutoff=0.2, 
                           status="hyper", dname = "sigDMR", geneset_cutoff=200, 
                           genome="hg38"){
    
  res <- rGREAT::great(gr, 
                       gene_sets = geneset, 
                       tss_source = paste0("GREAT:", genome),
                       mode = "oneClosest")
  
  genesetName <- gsub(":", "_", genesetName, fixed = TRUE)
  
  pdf(paste0("GREAT/distance_to_TSS_", status, "_", dname, ".pdf"))
  rGREAT::plotRegionGeneAssociations(res)
  rGREAT::plotVolcano(res)
  dev.off()
  
  
  res_table <- rGREAT::getEnrichmentTable(res) %>%
    dplyr::filter(p_adjust < padj_cutoff) |>
    dplyr::filter(gene_set_size < geneset_cutoff)
  
  genes <- sapply(res_table$id, function(x){
    gene_gr <- rGREAT::getRegionGeneAssociations(res, term_id = x, by_middle_points = FALSE,
                                         use_symbols = TRUE)
    gene <- unlist(gene_gr$annotated_genes)
    names(gene) <- NULL
    genes <- paste(unique(gene), collapse = ",")
  })
  res_table$gene_name <- genes
  
  write.table(res_table, 
              paste0("GREAT/enrichment_in_", genesetName, "_", status, "_", dname, ".tsv"), 
              row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
  
  gene_gr <- rGREAT::getRegionGeneAssociations(res, term_id = NULL, by_middle_points = FALSE,
                                             use_symbols = TRUE)
  gene <- gene_gr$annotated_genes
  gene_vec <- sapply(gene, function(g){
    names(g) <- NULL
    g <- paste(g, collapse = ",")
  })
  gene_gr$annotated_genes <- gene_vec
  dist <- gene_gr$dist_to_TSS
  dist_vec <- sapply(dist, function(d){
    d <- paste(d, collapse = " ")
  })
  gene_gr$dist_to_TSS <- dist_vec
  
  write.table(GenomicPlot::gr2df(gene_gr), 
              paste0("GREAT/Annotated_", genesetName, "_", status, "_", dname, ".bed"),
              row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
  
}


#' output_DMR
#' @title Output DMRs 
#' @description Write DMR to bed files
#' @param DMR a list of two dataframes.
#' @param dir path for saving DMR bed files.
#' 
#' @export output_DMR
#' 
output_DMR <- function(DMR){
  dmr <- DMR$sigRegions
  background <- DMR$bgRegions
  
  dmr$name <- rownames(dmr)
  dmr$strand <- rep("*", nrow(dmr))
  dmr <- dmr[, c("chr", "start", "end", "name", "stat", "strand", "length", "nCG", "ratio", "status")]
  
  write.table(dmr, "DMRs/DMR.bed", row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
  
  background$name <- rownames(background)
  background$strand <- rep("*", nrow(background))
  background <- background[, c("chr", "start", "end", "name", "stat", "strand", "length", "nCG", "ratio", "status")]
  
  write.table(background, "DMRs/background.bed", row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
  
  hyper <- dmr |>
    dplyr::filter(status == "hyper")
  hypo <- dmr |>
    dplyr::filter(status == "hypo")
  
  write.table(hyper, "DMRs/DMR_hyper.bed", row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
  write.table(hypo, "DMRs/DMR_hypo.bed", row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
}

#' DSS_pairwise
#' @title DML and DMR calling with DSS method for two groups
#' @description Call DML using two group design, then get DMRS and background regions.
#' @param bss a bsseq object.
#' @param pval_cutoff Numeric cutoff value [from 0 to 1] for the pval of DML used for DMR detection.
#' @param minDifff cutoff value [from 0 to inf] for the minimum difference between mean methylation 
#' levels between group1 and group2 during DMR detection.
#' @param minSites Numeric for the minimum number of Cytosines for a DMR.
#' @param condition1 Character indicating the group1. 
#' @param condition2 Character indicating the group2.
#' @param cores number of cpus.
#' @param ratio_cutoff cutoff value [from 0 to inf] for the ratio areaStat/nSites used for DMR detection.
#' 
#' @return a list of DMRs and background regions
#' 
#' @importFrom parallel detectCores
#' @importFrom DSS DMLtest callDML callDMR
#' 
#' @export DSS_pairwise

DSS_pairwise <- function(bss, condition1, condition2, pval_cutoff, minDiff, 
                         minSites=3, cores=5, ratio_cutoff=2){
  message("[DSS_pairwise] starting ... condition1: ", 
          condition1, " condition2: ", condition2)
  print(pData(bss))
  
  aname <- paste0(condition2, "_vs_", condition1)
  samples1 <- as.data.frame(pData(bss)) |>
    dplyr::filter(group == condition1) |>
    rownames()
  samples2 <- as.data.frame(pData(bss)) |>
    dplyr::filter(group == condition2) |>
    rownames()
 
  message("[DSS_pairwise] DML test ..")
  
  # Load the parallel package to make detectCores available
  # This is needed for DSS::DMLtest which uses detectCores internally
  library(parallel)
  
  DML= DSS::DMLtest(bss, 
               group1=samples1,
               group2=samples2, 
               smoothing=TRUE,
               ncores=cores)
  
  
  #message("[DSS_pairwise] call DML ..")
  
  #DMLs = DSS::callDML(DML, p.threshold = pval_cutoff, delta = 0.1)
  #DMLs <- na.omit(DMLs)
  
  message("[DSS_pairwise] call DMR ..")
  
  dmrs <- DSS::callDMR(DML, 
                  p.threshold = pval_cutoff, 
                  delta = 0.1, 
                  minlen=50, 
                  minCG=minSites, 
                  dis.merge=100, 
                  pct.sig = 0.5) |>
    dplyr::filter(abs(diff.Methy) > minDiff) |>
    dplyr::filter(!is.na(chr)) |>
    dplyr::mutate(ratio = areaStat/nCG) |>
    dplyr::filter(abs(ratio) > ratio_cutoff) |>
    dplyr::mutate(stat = areaStat, .keep = "unused") |>
    dplyr::mutate(status = case_when(stat > 0 ~ "hyper",
                                     stat < 0 ~ "hypo",
                                     .default = "none"))
  
  dmrs_background <- callDMR(DML, 
                  p.threshold = 1, 
                  delta = 0, 
                  minlen=50, 
                  minCG=minSites, 
                  dis.merge=0, # do not merge regions
                  pct.sig = 0.5) |> # keep same as dmrs to avoid merging regions
    dplyr::filter(!is.na(chr)) |>
    dplyr::mutate(ratio = areaStat/nCG) |>
    dplyr::mutate(stat = areaStat, .keep = "unused") |>
    dplyr::mutate(status = case_when(stat > 0 ~ "hyper",
                                     stat < 0 ~ "hypo",
                                     .default = "none"))
    
  
  message("[DSS_pairwise] finished ..")
  
  if(nrow(dmrs) > 0){
    return(list(sigRegions = dmrs, bgRegions = dmrs_background))
  }else{
    message("!!! No significant DMR found at p.threshold ", pval_cutoff)
    return(NULL)
  }
}
