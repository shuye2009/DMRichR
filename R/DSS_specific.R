
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
    bpparams <- BiocParallel::SnowParam(workers = detectCores()-2, progressbar = TRUE)
    smooth_params <- BiocParallel::SnowParam(workers = 1, progressbar = TRUE)
  }else{
    bpparams <- smooth_params <- BiocParallel::MulticoreParam(workers = cores, progressbar = FALSE)
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
  
  
  if(!file.exists("bss.RDS")){
    bss <- BSmooth(bs, BPPARAM = smooth_params)
    
    saveRDS(bss, "bss.RDS")
  }else{
    bss <- readRDS("bss.RDS")
  }
  
  return(bss)
}


#' DSS_multi_factor
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
#' @export DSS_multi_factor

DSS_multi_factor <- function(bss, design, factor1, factor2, pval_cutoff, ratio_cutoff, minSites=3){
  if(!file.exists("DMLfit.RDS")){
    if(is.null(factor2)){
      fomu <- as.formula(paste("~", factor1))
    }else{
      fomu <- as.formula(paste("~", factor1, "+", factor2, "+", factor1, ":", factor2))
    }
    DMLfit= DSS::DMLfit.multiFactor(bss, design = design, formula = fomu)
    saveRDS(DMLfit, "DMLfit.RDS")
  }else{
    DMLfit <- readRDS("DMLfit.RDS")
  }
  
  message("desing matrix:")
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
  gc()
  
  message("[DSS_multi_factor] finding DMR for factor ", factor1)
  
  dmrsfactor1 <- findDMR(DMLfactor1, 
                        pval_cutoff = pval_cutoff, 
                        ratio_cutoff = ratio_cutoff,
                        minSites = minSites)
  if(!is.null(factor2)){
    message("[DSS_multi_factor] finding DMR for factor ", factor2)
    dmrsfactor2 <- findDMR(DMLfactor2, 
                          pval_cutoff = pval_cutoff, 
                          ratio_cutoff = ratio_cutoff)
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
    dplyr::mutate(ratio = base::abs(areaStat/nCG)) |>
    dplyr::mutate(stat = areaStat, .keep = "unused") |>
    dplyr::mutate(status = case_when(stat > 0 ~ "hyper",
                                     stat < 0 ~ "hypo",
                                     .default = "none"))
  
  message("[findDMR] sample of significant DMR at pval_cutoff ", pval_cutoff, 
          " before filtering for ratio_cutoff ", ratio_cutoff)
  print(head(dmrs))
  
  if(nrow(dmrs) > 0){
    message("select significant DMR")
    dmrs <- dmrs |>
      dplyr::filter(ratio > ratio_cutoff)
  
    message("select background DMR")
    dmrs_background <- DSS::callDMR(DML, 
                               delta=0, 
                               p.threshold=0.99,
                               minlen=50, 
                               minCG=3, 
                               dis.merge=100, 
                               pct.sig=0.01) |>
      dplyr::filter(!is.na(chr)) |>
      dplyr::mutate(ratio = base::abs(areaStat/nCG)) |>
      dplyr::mutate(stat = areaStat, .keep = "unused") |>
      dplyr::mutate(status = case_when(stat > 0 ~ "hyper",
                                       stat < 0 ~ "hypo",
                                       .default = "none")) |>
      dplyr::filter(ratio < ratio_cutoff) |>
      rbind(dmrs)
    
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
#' @param padj_cutoff Numeric cutoff value [from 0 to 1] for the adjusted p-value.
#' @param geneset_cutoff cutoff value for the maximal size of the geneset.
#' @param status indicates the direction of methylation changes, either hyper or hypo.
#' @param dname indicates nature of DMR, either sigDMR or background
#' @param dir path for saving GREAT analysis results.
#' @param genome Character specifying the genome.
#' 
#' @export GREAT_analysis

GREAT_analysis <- function(gr, genesetName="GO:BP", padj_cutoff=0.2, 
                           status="hyper", dname = "sigDMR", geneset_cutoff=200, 
                           dir = "./", genome="hg38"){
    
  res <- rGREAT::great(gr, 
                       gene_sets = genesetName, 
                       tss_source = paste0("GREAT:", genome),
                       mode = "oneClosest")
  
  genesetName <- gsub(":", "_", genesetName, fixed = TRUE)
  
  pdf(file.path(dir, paste0("GREAT/distance_to_TSS_", status, "_", dname, ".pdf")))
  plotRegionGeneAssociations(res)
  plotVolcano(res)
  dev.off()
  
  
  res_table <- getEnrichmentTable(res) %>%
    dplyr::filter(p_adjust < padj_cutoff) |>
    dplyr::filter(gene_set_size < geneset_cutoff)
  
  genes <- sapply(res_table$id, function(x){
    gene_gr <- getRegionGeneAssociations(res, term_id = x, by_middle_points = FALSE,
                                         use_symbols = TRUE)
    gene <- unlist(gene_gr$annotated_genes)
    names(gene) <- NULL
    genes <- paste(unique(gene), collapse = ",")
  })
  res_table$gene_name <- genes
  
  write.table(res_table, 
              file.path(dir, paste0("GREAT/enrichment_in_", genesetName, "_", status, "_", dname, ".tsv")), 
              row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
  
  gene_gr <- getRegionGeneAssociations(res, term_id = NULL, by_middle_points = FALSE,
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
              file.path(dir, paste0("GREAT/Annotated_", genesetName, "_", status, "_", dname, ".bed")),
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
output_DMR <- function(DMR, dir){
  dmr <- DMR$sigRegions
  background <- DMR$bgRegions
  
  dmr$name <- rownames(dmr)
  dmr$strand <- rep("*", nrow(dmr))
  dmr <- dmr[, c("chr", "start", "end", "name", "stat", "strand", "length", "nCG", "status")]
  
  write.table(dmr, file.path(dir,"DMR/DMR.bed"), row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
  
  background$name <- rownames(background)
  background$strand <- rep("*", nrow(background))
  background <- background[, c("chr", "start", "end", "name", "stat", "strand", "length", "nCG", "status")]
  
  write.table(background, file.path(dir, "DMR/background.bed"), row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
  
  hyper <- dmr |>
    dplyr::filter(status == "hyper")
  hypo <- dmr |>
    dplyr::filter(status == "hypo")
  
  write.table(hyper, file.path(dir,"DMRs/DMR_hyper.bed"), row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
  write.table(hypo, file.path(dir,"DMRs/DMR_hypo.bed"), row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
}
