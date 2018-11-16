# DMRichR
#### A workflow for the statistical analysis and visualization of differentially methylated regions (DMRs) of CpG count matrices (Bismark cytosine reports) from the [CpG_Me pipeline](https://github.com/ben-laufer/CpG_Me).

### Table of Contents
1. [Installation](https://github.com/ben-laufer/DMRichR#installation)
2. [The Design Matrix and Covariates](https://github.com/ben-laufer/DMRichR#the-design-matrix-and-covariates)
3. [Input](https://github.com/ben-laufer/DMRichR#input)
   1. [Generic Example](https://github.com/ben-laufer/DMRichR#generic-example)
   2. [UC Davis Example](https://github.com/ben-laufer/DMRichR#uc-davis-example)
4. [Output](https://github.com/ben-laufer/DMRichR#output)
5. [DMR Interpretation](https://github.com/ben-laufer/DMRichR#dmr-interpretation)
6. [Citation](https://github.com/ben-laufer/DMRichR#citation)
7. [Acknowledgements](https://github.com/ben-laufer/DMRichR#acknowledgements)

## Installation

No manual installation of R packages is required, since the required packages and updates will occur automatically upon running the [executable script](exec/DM.R) located in the `exec` folder. However, the package does require Bioconductor 3.8, which you can install or update to using:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("BiocInstaller", version = "3.8")
```

Additionally, if you are interested in creating your own workflow as opposed to using the executable script, you can download the package using:

`BiocManager::install(c("remotes", "ben-laufer/DMRichR"))`

## The Design Matrix and Covariates

This script requires a basic design matrix to identify the groups and covariates, which should be named `sample_info.csv` and contain header columns to identify the factor. It is important to have the label for the experimental samples start with a letter in the alphabet that comes after the one used for control samples in order to obtain results for experimental vs. control rather than control vs. experimental. You can select which specific samples to analzye from the working directory through the design matrix, where pattern matching of the sample name will only select bismark cytosine report files with a matching name before the first underscore. Within the script, covariates can be selected for adjustment. There are two different ways to adjust for covariates: directly adjust values or balance permutations.


| Name          | Diagnosis      | Age           |  Sex          |
| ------------- | -------------- | ------------- | ------------- |
| SRR3537014    | Idiopathic_ASD | 14            | M             |
| SRR3536981    | Control        | 42            | F             |


## Input

Before running the executable, ensure you have the following project directory tree structure for the [Bismark cytosine reports](https://github.com/ben-laufer/CpG_Me) and design matrix:

```
├── Project
│   ├── cytosine_reports
│   │   ├── sample1_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz
│   │   ├── sample2_bismark_bt2.deduplicated.bismark.cov.gz.CpG_report.txt.gz
│   │   ├── sample_info.csv
```

This workflow requires the following variables:
1. `-g --genome` Select either: hg38, mm10, rn6, or rheMac8
2. `-x --coverage` The coverage cutoff for all samples, 1x is recommended
3. `-t --testCovariate` The covariate to test for significant differences between experimental and control, i.e.: Diagnosis
4. `-a --adjustCovariate` Adjust covariates that are continuous or contain two or more groups. More than one covariate can be adjusted for., i.e.: "Age" or c("Age", "PMI")
5. `-m --matchCovariate` Covariate to balance permutations, which is ideal for two group covariates. Only one covariate can be balanced. i.e: Sex
6. `-c --cores` The number of cores to use, 1 is recommended during this stage of development.

#### Generic Example

Below is an example of how to execute the [main R script (DM.R)](exec/DM.R) in the `exec` folder on command line. This should be called from the working directory that contains the cytosine reports.

```
call="Rscript \
--vanilla \
/share/lasallelab/programs/DMRichR/DM.R \
--genome hg38 \
--coverage 1 \
--testCovariate Diagnosis \
--adjustCovariate Age \
--matchCovariate Sex \
--cores 1"

echo $call
eval $call
```

#### UC Davis Example

If you are using the Barbera cluster at UC Davis, the following commands can be used before the above to execute `DM.R` from your login node (i.e. epigenerate), where `htop` should be called first to make sure the resources (2 cores and 64 GB RAM for a few days) are available. This should be called from the working directory that contains the cytosine reports and from within a `screen`, where you can use <kbd>⌃ Control</kbd> + <kbd>a</kbd> followed by <kbd>⌃ Control</kbd> + <kbd>d</kbd> to deatch after launching and reattach via `screen -r DMRs`.

```
screen -s DMRs
kinit -l 10d
module load R
```

## Output

This workflow provides the following files:
1. CpG methylation and coverage value distribution plots
2. DMRs and background regions
3. Individual smoothed methylation values for DMRs, background regions, and windows/bins
4. Smoothed global and chromosomal methylation values and statistics
5. Heatmap of DMRs
6. PCA plots of 20 Kb windows (all genomes) and CpG island windows (hg38, mm10, and rn6)
7. Gene ontology and pathway enrichments (enrichr for all genomes, GREAT for hg38 and mm10)
8. Gene region and CpG annotations and plots (hg38, mm10, or rn6)
9. Manhattan and Q-Qplots 
10. Blocks of methylation and background blocks

## DMR Interpretation

Fold changes are not utilized in this workflow. Rather, the focus is the beta coefficient, which is representative of the average [effect size](https://www.leeds.ac.uk/educol/documents/00002182.htm); however, it is on the scale of the [arcsine transformed differences](https://www.ncbi.nlm.nih.gov/pubmed/29481604) and must be divided by π (3.14) to be similar to the mean methylation difference over a DMR, which is provided in the `percentDifference` column. There ~~is~~ *was* also the raw difference column, which shows the percent difference in raw (non-smoothed methylation values) and in general will reflect overall difference in a DMR, although in some cases the weighting of the smoothing may change the value and its directionality. You can also read a general summary of the drmseq approach on [EpiGenie](https://epigenie.com/dmrseq-powers-whole-genome-bisulfite-sequencing-analysis/).

## Citation

If you use **DMRichR** in published research please cite the following 3 articles:

Laufer BI, Hwang H, Vogel Ciernia A, Mordaunt CE, LaSalle JM. Whole genome bisulfite sequencing of Down syndrome brain reveals regional DNA hypermethylation and novel disease insights. *bioRxiv*, 2018. **doi**: [10.1101/428482](https://doi.org/10.1101/428482)

Korthauer K, Chakraborty S, Benjamini Y, and Irizarry RA. Detection and accurate false discovery rate control of differentially methylated regions from whole genome bisulfite sequencing. *Biostatistics*, 2018. **doi**: [10.1093/biostatistics/kxy007](https://doi.org/10.1093/biostatistics/kxy007)

Hansen KD, Langmead B, Irizarry RA. BSmooth: from whole genome bisulfite sequencing reads to differentially methylated regions. *Genome Biology*, 2012. **doi**: [10.1186/gb-2012-13-10-r83](https://doi.org/10.1186/gb-2012-13-10-r83)

## Acknowledgements

This workflow is primarily based on the [dmrseq](https://www.bioconductor.org/packages/release/bioc/html/dmrseq.html) and [bsseq](https://www.bioconductor.org/packages/release/bioc/html/bsseq.html) bioconductor packages. I would like to thank [Keegan Korthauer](https://github.com/kdkorthauer), the creator of dmrseq, for helpful conceptual advice in establishing and optimizing this workflow. I would like to thank [Matt Settles](https://github.com/msettles) from the [UC Davis Bioinformatics Core](https://github.com/ucdavis-bioinformatics) for advice on creating an R package and use of the tidyverse. I would like to thank Rochelle Coulson for a script that was developed into the PCA function. I would also like to thank Blythe Durbin-Johnson and Annie Vogel Ciernia for statistical consulting that enabled the global and chromosomal methylation statistics. Finally, I would like to thank [Nikhil Joshi](https://github.com/najoshi) from the [UC Davis Bioinformatics Core](https://github.com/ucdavis-bioinformatics) for troubleshooting of a [resource issue](https://github.com/kdkorthauer/dmrseq/commit/38dea275bb53fcff3a0df93895af759b15c90e3e).
