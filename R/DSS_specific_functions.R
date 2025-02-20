

processReport <- function(report_path, file_pattern, design, cores){
  files <- list.files(report_path, file_pattern)
  if(Sys.info()['sysname'] == "Windows"){
    bpparams <- BiocParallel::SnowParam(workers = detectCores()-2, progressbar = TRUE)
    smooth_params <- BiocParallel::SnowParam(workers = 1, progressbar = TRUE)
  }else{
    bpparams <- smooth_params <- BiocParallel::MulticoreParam(workers = cores, progressbar = FALSE)
  }
  
  print(glue::glue("Reading cytosine reports..."))
  if(!file.exists("bs.RDS")){
    bs <- bsseq::read.bismark(files = file.path(report_path, files),
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

DSS_multi_factor <- function(bss, desgin, factor1, factor2, pval_cutoff, ratio_cutoff){
  if(!file.exists("DMLfit.RDS")){
    fomu <- as.formula(paste("~", factor1, "+", factor2, "+", factor1, ":", factor2))
    DMLfit= DSS::DMLfit.multiFactor(bss, design = design, formula = fomu)
    saveRDS(DMLfit, "DMLfit.RDS")
  }else{
    DMLfit <- readRDS("DMLfit.RDS")
  }
  
  DMLfit$X                 
  
  DMLfactor1 <- DSS::DMLtest.multiFactor(DMLfit, coef = 2)
  DMLfactor2 <- DSS::DMLtest.multiFactor(DMLfit, coef = 3)
  DMLInter <- DSS::DMLtest.multiFactor(DMLfit, coef = 4)
  
  rm(DMLfit)
  gc()
  
  
  pval_cutoff <- 0.05
  ratio_cutoff <- 2
  
  dmrsfactor1 <- findDMR(DMLfactor1, 
                                  pval_cutoff = pval_cutoff, 
                                  ratio_cutoff = ratio_cutoff)
  dmrsfactor2 <- findDMR(DMLfactor2, 
                                  pval_cutoff = pval_cutoff, 
                                  ratio_cutoff = ratio_cutoff)
  dmrsInter <- findDMR(DMLInter, 
                                  pval_cutoff = pval_cutoff, 
                                  ratio_cutoff = ratio_cutoff)
  return(list("factor1"= dmrsfactor1, "factor2"= dmrsfactor2, "interaction" = dmrsInter))
}

findDMR <- function(DML, pval_cutoff=0.05, ratio_cutoff=2){
  dmrs <- DSS::callDMR(DML, 
                  delta=0, 
                  p.threshold=pval_cutoff,
                  minlen=50, 
                  minCG=3, 
                  dis.merge=100, 
                  pct.sig=0.5) |>
    dplyr::mutate(stat = areaStat) |>
    dplyr::filter(abs(stat/nCG) > ratio_cutoff) |>
    dplyr::mutate(status = case_when(stat > 0 ~ "hyper",
                                     stat < 0 ~ "hypo",
                                     .default = "none"))
  
  dmrs_background <- DSS::callDMR(DML, 
                             delta=0, 
                             p.threshold=1.0,
                             minlen=50, 
                             minCG=3, 
                             dis.merge=100, 
                             pct.sig=0.5) |>
    dplyr::mutate(stat = areaStat) |>
    dplyr::filter(abs(stat/nCG) < ratio_cutoff) |>
    dplyr::filter(length < max(dmrs$length)) |>
    dplyr::mutate(status = case_when(stat > 0 ~ "hyper",
                                     stat < 0 ~ "hypo",
                                     .default = "none")) |>
    rbind(dmrs)
  
  
  
  return(list(sigRegions=dmrs, bgRegions=dmrs_background))
}

GREAT_analysis <- function(DMR, genesetName="GO:BP", padj_cutoff=0.2, geneset_cutoff=200){
  dir.create(("GREAT"))
  
  dmr_list <- list()
  dmr_list$sigDMR <- DMR$sigRegions
  dmr_list$background <- dmrs$bgRegions
  
  for(dname in names(dmr_list)[1]){
    dmr <- dmr_list[[dname]]
    
    hyper_gr <- as(dmr[dmr$status=="hyper",], "GRanges")
    hypo_gr <- as(dmr[dmr$status=="hypo",], "GRanges")
    
    #genesetName <- "GO:BP"
    hyper_res <- rGREAT::great(hyper_gr, gene_sets = genesetName, "GREAT:hg38")
    hypo_res <- rGREAT::great(hypo_gr, gene_sets = genesetName, "GREAT:hg38")
    
    genesetName <- gsub(":", "_", genesetName, fixed = TRUE)
    
    pdf(paste0("GREAT/distance_to_TSS_hyper_", dname, ".pdf"))
    plotRegionGeneAssociations(hyper_res)
    plotVolcano(hyper_res)
    dev.off()
    
    pdf(paste0("GREAT/distance_to_TSS_hypo_", dname, ".pdf"))
    plotRegionGeneAssociations(hypo_res)
    plotVolcano(hypo_res)
    dev.off()
    
    hyper_res_table <- getEnrichmentTable(hyper_res) %>%
      dplyr::filter(p_adjust < padj_cutoff) |>
      dplyr::filter(gene_set_size < geneset_cutoff)
    
    hypo_res_table <- getEnrichmentTable(hypo_res) %>%
      dplyr::filter(p_adjust < padj_cutoff) |>
      dplyr::filter(gene_set_size < geneset_cutoff)
    
    hyper_genes <- sapply(hyper_res_table$id, function(x){
      gene_gr <- getRegionGeneAssociations(hyper_res, term_id = x, by_middle_points = FALSE,
                                           use_symbols = TRUE)
      gene <- unlist(gene_gr$annotated_genes)
      names(gene) <- NULL
      genes <- paste(unique(gene), collapse = ",")
    })
    hyper_res_table$gene_name <- hyper_genes
    
    hypo_genes <- sapply(hypo_res_table$id, function(x){
      gene_gr <- getRegionGeneAssociations(hypo_res, term_id = x, by_middle_points = FALSE,
                                           use_symbols = TRUE)
      gene <- unlist(gene_gr$annotated_genes)
      names(gene) <- NULL
      genes <- paste(unique(gene), collapse = ",")
    })
    hypo_res_table$gene_name <- hypo_genes
    
    write.table(hyper_res_table, 
                paste0("GREAT/enrichment_in_", genesetName, "_hyper_", dname, ".tsv"), 
                row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
    write.table(hypo_res_table, 
                paste0("GREAT/enrichment_in_", genesetName, "_hypo_", dname, ".tsv"), 
                row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
    
    hyper_gene_gr <- getRegionGeneAssociations(hyper_res, term_id = NULL, by_middle_points = FALSE,
                                               use_symbols = TRUE)
    hyper_gene <- hyper_gene_gr$annotated_genes
    hyper_gene_vec <- sapply(hyper_gene, function(g){
      names(g) <- NULL
      g <- paste(g, collapse = ",")
    })
    hyper_gene_gr$annotated_genes <- hyper_gene_vec
    hyper_dist <- hyper_gene_gr$dist_to_TSS
    hyper_dist_vec <- sapply(hyper_dist, function(d){
      d <- paste(d, collapse = " ")
    })
    hyper_gene_gr$dist_to_TSS <- hyper_dist_vec
    
    write.table(GenomicPlot::gr2df(hyper_gene_gr), 
                paste0("GREAT/Annotated_", genesetName, "_hyper_", dname, ".bed"),
                row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
    
    hypo_gene_gr <- getRegionGeneAssociations(hypo_res, term_id = NULL, by_middle_points = FALSE,
                                              use_symbols = TRUE)
    hypo_gene <- hypo_gene_gr$annotated_genes
    hypo_gene_vec <- sapply(hypo_gene, function(g){
      names(g) <- NULL
      g <- paste(g, collapse = ",")
    })
    hypo_gene_gr$annotated_genes <- hypo_gene_vec
    hypo_dist <- hypo_gene_gr$dist_to_TSS
    hypo_dist_vec <- sapply(hypo_dist, function(d){
      d <- paste(d, collapse = " ")
    })
    hypo_gene_gr$dist_to_TSS <- hypo_dist_vec
    
    write.table(GenomicPlot::gr2df(hypo_gene_gr), 
                paste0("GREAT/Annotated_", genesetName, "_hypo_", dname, ".bed"),
                row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
  }
  
}

output_DMR <- function(DMR, dir){
  dir.create(file.path(dir,"DMR"))
  dmr <- DMR$sigRegions
  background <- dmrs$bgRegions
  
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
  
  write.table(hyper, file.path(dir,"DMR/DMR_hyper.bed"), row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
  write.table(hypo, file.path(dir,"DMR/DMR_hypo.bed"), row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
}
