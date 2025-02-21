
create_design <- function(report_path, report_pattern, exclude = NA, outdir="."){
  files <- list.files(path = report_path, pattern = report_pattern)
  # remove batch-0
  unless(is.na(exclude)) files <- files[!grepl(exclude, files)]
  names <- gsub(report_pattern, "", files)
  
  design <- as.data.frame(Reduce(rbind, strsplit(names, split="-"))) |>
    mutate(name = files)
  colnames(design) <- c("pheno", "treat", "rep", "name")
  rownames(design) <- names
  
  design <- design |>
    mutate(treat=factor(treat, levels=c("NIR", "IR"))) |>
    mutate(pheno=factor(pheno, levels= c("WT", "R172K"))) |>
    mutate(group=factor(paste(pheno, treat, sep="-"))) |>
    mutate(group=relevel(group, "WT-NIR"))
  write.table(design, file.path(outdir, "sample_info.txt"), col.names = T, row.names = T, sep="\t")
}

outdir <- "/cluster/projects/hardinggroup/Shuye/EMseq/Harding_combined/report/DSS_analysis"
reportPath <- "/cluster/projects/hardinggroup/Shuye/EMseq/Harding_combined/report"
reportPattern <- ".CG_report.txt.gz"

create_design(report_path, report_pattern, exclude = "-0", outdir)