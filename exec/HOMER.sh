#!/bin/bash

# HOMER transcription factor motif analysis for DMRichR
# By Ben Laufer

homerDMR(){
	genome=${1}
	cores=${2}
	
	cd HOMER
	
    echo 
    echo "Testing both hypermethylated and hypomethylated DMRs"
	mkdir both
	
	call="findMotifsGenome.pl \
	DMRs.bed \
	${genome} \
	both/ \
	-bg background.bed \
	-cpg \
	-size given \
	-p ${cores} \
	-nomotif"

	echo $call
	eval $call
    
    echo 
	echo "Testing hypermethylated DMRs"
	mkdir hyper

	call="findMotifsGenome.pl \
	DMRs_hyper.bed \
	${genome} \
	hyper/ \
	-bg background.bed \
	-cpg \
	-size given \
	-p ${cores} \
	-nomotif"

	echo $call
	eval $call

    echo 
    echo "Testing hypomethylated DMRs"
	mkdir hypo

	call="findMotifsGenome.pl \
	DMRs_hypo.bed \
	${genome} \
	hypo/ \
	-bg background.bed \
	-cpg \
	-size given \
	-p ${cores} \
	-nomotif"

	echo $call
	eval $call
	
## annotate DMRs, Added by Shuye Pu

  echo 
      echo "Annotating DMRs"
  call="annotatePeaks.pl \
	DMRs.bed \
	${genome} \
	-annStats both/DMRs_annot_stat.tab \
	-p ${cores} \
	> both/DMRs_annot.tab"

	echo $call
	eval $call
    
  call="annotatePeaks.pl \
	DMRs_hyper.bed \
	${genome} \
	-annStats hyper/DMRs_hyper_annot_stat.tab \
	-p ${cores} \
	> hyper/DMRs_hyper_annot.tab"

	echo $call
	eval $call

 call="annotatePeaks.pl \
	DMRs_hypo.bed \
	${genome} \
	-annStats hypo/DMRs_hypo_annot_stat.tab \
	-p ${cores} \
	> hypo/DMRs_hypo_annot.tab"

	echo $call
	eval $call

}
export -f homerDMR

homerDMR ${1} ${2}
