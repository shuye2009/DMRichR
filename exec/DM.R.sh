#!/bin/bash
#
#SBATCH --job-name=DM.R
#SBATCH --ntasks=20 # Number of cores/threads
#SBATCH --mem=64000 # Ram in Mb
#SBATCH --partition=production 
#SBATCH --time=2-00:00:00

##########################################################################################
# Author: Ben Laufer
# Email: blaufer@ucdavis.edu 
##########################################################################################

###################
# Run Information #
###################

start=`date +%s`

hostname

THREADS=${SLURM_NTASKS}
MEM=$(expr ${SLURM_MEM_PER_CPU} / 1024)

echo "Allocated threads: " $THREADS
echo "Allocated memory: " $MEM

################
# Load Modules #
################

module load R
module load homer

########
# DM.R #
########

call="Rscript \
--vanilla \
/share/lasallelab/programs/DMRichR/DM.R \
--genome hg38 \
--coverage 1 \
--perGroup '0.75' \
--minCpGs 5 \
--maxPerms 10 \
--maxBlockPerms 10 \
--cutoff '0.05' \
--testCovariate Diagnosis \
--adjustCovariate 'Sex;Age' \
--sexCheck TRUE \
--GOfuncR TRUE \
--EnsDb FALSE \
--cores 20 \
--filePattern '*.CG_report.txt.gz' \
--internet FALSE"

echo $call
eval $call

###################
# Run Information #
###################

end=`date +%s`
runtime=$((end-start))
echo $runtime
