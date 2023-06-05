#!/bin/bash

# batch_qc_sfs.sh

BAMPREFIX='/space/s1/linderoth/imitator/spades/imitator'
SITESPREFIX='imi_model_allsites'
OUTPREFIX='/space/s2/diana/Rimitator/sites_290421/sfs/'
REF='/space/s1/linderoth/imitator/spades/ref/imi_combined_targetedAndflanking_geneid.fasta'
ANGSD='/home/tyler/install/angsd'
MORPHS=(banded striped admixed)
SITES=(set2 set1 set3)
declare -a SAFIDX=()
declare -a SFSOUT=()
declare -a LOGS=()

for pos in "${SITES[@]}"
do
	POSLIST="${SITESPREFIX}_${pos}.pos"
	RFLIST="${SITESPREFIX}_${pos}.rf"
	for pop in "${MORPHS[@]}"
	do
		BAM="${BAMPREFIX}_${pop}_bams.txt"
		OUT="${OUTPREFIX}${pop}_allsites_${pos}"
		OUTLOG="${OUTPREFIX}${pop}_allsites_${pos}_sfs.log"
		"${ANGSD}/angsd" -bam $BAM -sites $POSLIST -rf $RFLIST -GL 1 -doSaf 1 -fold 1 -minQ 20 -minMapQ 20 -anc $REF -out $OUT -P 4 2> $OUTLOG &
		SAFIDX+=("${OUT}.saf.idx")
		SFSOUT+=("${OUT}.sfs")
		LOGS+=($OUTLOG)
	done
	wait
done
wait

i=0
for idx in "${SAFIDX[@]}"
do
	"${ANGSD}/misc/realSFS" $idx > ${SFSOUT[$i]} 2> ${LOGS[$i]} &
	((i++))
done
wait

exit
