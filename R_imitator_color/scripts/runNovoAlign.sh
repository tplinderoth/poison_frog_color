#!/bin/bash

# runNovoAlign.sh
# use samtools 1.3-20-gd49c73b (using htslib 1.3-34-g820689c)
# use novoalign V 3.04.06

### define variables ###
VERSION=0.0.6
FQDIR="" # directory containing fastq files
REF="" # file containing reference fasta -d
OUTDIR="" # output directory
INDEXNAME="" # reference index name
FQFMT='STDFQ' # format of fastq files (quality score offset) -F
MAXSCORE="" # alignment score threshold -t
READGROUP="" # read group for reads
#######################################################################
# use default alignment score threshold, -t log4(N),4.5
# -t A,B
# threshold = (L-A)*B, L=read length
# without specifying threshold, novoalign will calculate a threshold
# for each read such that P(finding false positive alignment) < 0.001.
# A mismatch at a base with high phred quality will score 30 points so
# a threshold of 90 would allow at least 3 mismatches.
# alingment score is -10log10(P(R|Ai)), where P(R|Ai) is prob of read
# sequence given alignment location i.
#######################################################################
#GAPOPEN # gap opening penalty -t
GAPEXTEND=6 # gap extend penalty -x (novoaling default = 6)
MATCHRWD=$((GAPEXTEND/2)) # sets a match reward, suggested is ~half gap extend penalty, --matchreward
MINBASE=0 # minimum number of good quality base pairs for a read -l, should be ~half read length
#ADAPTER1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' # adapter sequence 1 -a # only available in licensed version
#ADAPTER2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA' # adapter sequence 2 -a # only avilable in licensed version
MAPDIFF=5 # score difference between first two alignments for reporting repeats -R
REPEAT='None' # method for handling reads with multiple mapping locations -r
REPLIMIT=10000 # limit for the number of repeat alignments reported for Exhaustive method
AVGSZ=250 # average insert size
SDSZ=50 # standard deviation for insert size
MINMAPQ=20 # minimum mapping quality to retain read -10*log10(P(mapping position is wrong))

if [ "$#" -lt 10 ]
then
>&2 echo -e "\nrunNovoAlign V $VERSION"
>&2 echo -e "Wrapper around NovoAlign and SAMtools to produce BAM files from paired-end fastq files"
>&2 echo -e "\nDependencies:\nNovoAlign\nSAMtools"
>&2 echo -e "\nUsage: bash runNovoAlign.sh [options]"
>&2 echo -e "\nOPTIONS:"
>&2 echo -e "-f|--fqdir       <string>   Directory with fastq files [$FQDIR]"
>&2 echo -e "-r|--ref         <string>   Reference fasta file [$REF]"
>&2 echo -e "-o|--outdir      <string>   Output directory [$OUTDIR]"
>&2 echo -e "-n|--idxname     <string>   Name of indexed reference (optional) [$INDEXNAME]"
>&2 echo -e "-g|--group       <string>   Read group name for SAM @RG ID flag [$READGROUP]"
>&2 echo -e "-b|--minbase     <int>      Minimum number of good quality read bases [$MINBASE]"
>&2 echo -e "-s|--maxscore    <int>      Highest acceptable alignment score [novoalign default]"
>&2 echo -e "-d|--scorediff   <int>      Max alignment score difference for reporting repeats [$MAPDIFF]"
>&2 echo -e "-t|--repmethod   <string>   Method for handling repeats: 'None', 'Random', 'All', 'Exhaustive' [$REPEAT]"
>&2 echo -e "-m|--minmapq     <int>      Minimum mapping quality for retaining reads [$MINMAPQ]"
>&2 echo -e "-s|--avgsize     <int>      Average insert size [$AVGSZ]"
>&2 echo -e "-v|--stdev       <int>      Insert size standard deviation [$SDSZ]"
>&2 echo -e ""
exit
fi

### parse user arguments ###

while [[ $# -gt 1 ]]
	do
	key="$1"

case $key in
	-f|--fqdir)
	FQDIR="$2"
	shift
	;;
	-r|--ref)
	REF="$2"
	shift
	;;
	-o|--outdir)
	OUTDIR="$2"
	shift
	;;
	-n|--idxname)
	INDEXNAME="$2"
	shift
	;;
	-g|--group)
	READGROUP="$2"
	shift
	;;
	-b|--minbase)
	MINBASE="$2"
	shift
	;;
	-s|--maxscore)
	MAXSCORE="$2"
	shift
	;;
	-d|--scorediff)
	MAPDIFF="$2"
	shift
	;;
	-t|--repmethod)
	REPEAT="$2"
	shift
	;;
	-m|--minmapq)
	MINMAPQ="$2"
	shift
	;;
	-s|--avgsize)
	AVGSZ="$2"
	shift
	;;
	-v|--stdev)
	SDSZ="$2"
	shift
	;;
	*)
	>&2 echo "Unknown command $1"
	exit 1
	;;
esac
shift
done

### perform input checks ###

>&2 echo

# check fastq directory
if [ -d "$FQDIR" ]
then
	if [[ $FQDIR == */ ]]
	then
		FQDIR=$( echo "$FQDIR" | sed 's/\/$//' )
	fi
else
	>&2 echo "Error: directory containing fastq files does not exist --> exiting"
	exit 1;
fi
>&2 echo "Fastq directory: $FQDIR"

# check indexed reference
if [ -f "$REF" ]
then
	if [ ! -s "$REF" ]
	then
		>&2 echo "Error: reference file has zero size --> exiting"
		exit 1
	fi
else
	>&2 echo "Error: index file does not exist --> exiting\n"
	exit 1
fi
>&2 echo "Reference: $REF"
refdir=$( dirname "${REF}"  )

# check for indexed reference name
if [ -z "$INDEXNAME" ]
then
	INDEXNAME=${REF%.*}; INDEXNAME=${INDEXNAME#/}; INDEXNAME=$( echo "$INDEXNAME" | sed 's/[^/\]*\///g' )
else
	if [[ $INDEXNAME == /* ]]
	then
		INDEXNAME=${INDEXNAME#/}; INDEXNAME=$( echo "$INDEXNAME" | sed 's/[^/\]*\///g' )
	fi
fi
INDEXNAME="${refdir}/${INDEXNAME}"
>&2 echo -e "Indexed Ref Name: $INDEXNAME"

# make sure index name is different than reference to prevent overwrite
refbase=$( basename "${REF}" )
if [ "$INDEXNAME" = "${refdir}/${refbase}"  ]
then
	>&2 echo -e "\nreference index name must be different than reference name --> exiting\n"
	exit 1
fi

# check output directory
if [ -d "$OUTDIR" ]
then
	if [[ $OUTDIR == */ ]]
	then
		OUTDIR=$( echo "$OUTDIR" | sed 's/\/$//' )
	fi
else
	>&2 echo "Error: output directory does not exist --> exiting"
	exit 1
fi
>&2 echo -e "Dumping BAMS to: $OUTDIR\n"

# check fastq format

case "$FQFMT" in
STDFQ)
	FMTMESSAGE="Sanger coding quality scores"
	;;
SLXFQ)
	FMTMESSAGE="Solexa style quality scores"
	;;
ILMFQ)
	FMTMESSAGE="Illumina coding quality scores"
	;;
*)
	FMTMESSAGE="an invalid quality score encoding --> exiting"
	exit 1
	;;
esac

>&2 echo "Fastq file format set to $FMTMESSAGE"

# check gap extend penalty
if [ "$GAPEXTEND" -lt 0 ]
then
	>&2 echo "Error: Gap extend penalty must be >= 0"
	exit 1
fi

# check match reward
if [ "$MATCHRWD" -lt 0 ]
then
	>&2 echo "Error: Match reward must be >= 0"
	exit 1
fi

# check max alignment score

#if [ -n "$MAXSCORE" ]
#then
#	if [ "$MAXSCORE" -le 0 ]
#	then
#		>&2 echo "Error: Maximum alignment score must be positive --> exiting"
#		exit 1
#	else
#		>&2 echo "Maximum alignment score: $MAXSCORE"
#	fi
#else
#	>&2 echo "Novoalign will use default score threshold"
#fi

# check minimum number of good quality read bases for mappping

if [ -z "$MINBASE" ]
then
	>&2 echo "Error: Must specify minimum number of good quality bases for read with -b || --minbase"
else
	if [ "$MINBASE" -lt 0 ]
	then
		>&2 echo "Error: Minimum number of good quality read bases must be >= 0"
	else
		>&2 echo "Minimum number of good quality read bases: $MINBASE"
	fi
fi

# check multiple mapping location cutoff

if [ "$MAPDIFF" -le 0 ]
then
	>&2 echo "Error: Score for indicating repeats must be positive --> exiting"
	exit 1
else
	>&2 echo "Repeat score difference threshold: $MAPDIFF"
fi

# check method for handling repeats
if [ $REPEAT == 'None' ]
then
	>&2 echo "No repeat alignments will be reported"
elif [ $REPEAT == 'Random' ]
then
	>&2 echo "A single random repeat alignment will be reported"
elif [ $REPEAT == 'All' ]
then
	>&2 echo "All repeats with score less than best alignment + $MAPDIFF will be reported"
elif [ $REPEAT == 'Exhaustive' ]
then
	REPEAT="E $REPLIMIT"
	>&2 echo "All repeats with P(R|Ai) score <= $MAPDIFF + $MAPDIFF will be reported"
else
	>&2 echo "$REPEAT is an invalid method for treating repeats: check -t/--repmethod"
	exit 1
fi

# check minimum mapping quality

if [ "$MINMAPQ" -lt 0 ]
then
	>&2 echo "Error: Minimum mapping quality must be >= 0 --> exiting"
	exit 1
else
	>&2 echo "Minimum Mapping Quality: $MINMAPQ"
fi

# check insert size and standard deviation

if [ "$AVGSZ" -le 0 ]
then
	>&2 echo "Error: Average insert size must be positive --> exiting"
	exit 1
	if [ "$SDSZ" -lt 0 ]
	then
		>&2 echo "Error: Standard deviation for insert size less than zero --> exiting"
		exit 1
	fi
else
	>&2 echo -e "\nInsert:\nAverage Size = $AVGSZ bp\nStandard Deviation = $SDSZ bp\n"
fi

# print adapter information
#>&2 echo -e "Trim 3' adapters:\nadapter1: $ADAPTER1\nadapter2: $ADAPTER2\n"

### index reference if not already done ###

if [ ! -e "$INDEXNAME" ]
then
	novoindex $INDEXNAME $REF
	rc=$?; if [[ $rc != 0 ]]; then echo "Call to Novoindex failed with return code $rc --> exiting"; exit $rc; fi
fi
if [ -e "$INDEXNAME" ] && [ ! -s "$INDEXNAME" ]
then
	>&2 echo "Indexed reference file $INDEXNAME appears to be empty"
	exit 1
fi

### process fastq files ###

for read1 in ${FQDIR}/*_R1.fastq
do

# collect file name information

read2=$( echo "$read1" | sed 's/_R1.fastq/_R2.fastq/' )
readu=$( echo "$read1" | sed 's/_R1.fastq/_u.fastq/' )

fullname=${read1%_clean_R1*}; fullname=${fullname#/}; fullname=$( echo "$fullname" | sed 's/[^\/]*\///g' )
if [ -z "$fullname" ]
then
	>&2 echo "Failed to extract fastq name info from $read1"
	exit 1
fi

lib=$( echo "$fullname" | sed 's/\([^_]\+\)_.*/\1/' )
if [ -z "$lib" ]
then
	>&2 echo "Couldn't determine library name for $read1"
	exit 1
fi

sample=$( echo "$fullname" | sed 's/[^_]\+_\([^_]\+\)_*.*/\1/' )
if [ -z "$sample" ]
then
	>&2 echo "Couldn't determine sample name for $read1"
	exit 1
fi

>&2 echo -e "\nProcessing library $fullname\n"

# run novoalign

if [ -z "$MAXSCORE" ]
then
>&2 echo "Novoalign will use default score threshold"

novoalign -o SAM "@RG\tID:${READGROUP}\tSM:${sample}\tLB:${lib}" -d "$INDEXNAME" -f "$read1" "$read2" -F "$FQFMT" -R "$MAPDIFF" -r "$REPEAT" -i PE "$AVGSZ","$SDSZ" -x "$GAPEXTEND" -l "$MINBASE" > "${OUTDIR}/NovoPaired.sam"
rc=$?; if [[ $rc != 0 ]]; then echo "Call to Novoalign failed on paired library $fullname with return code $rc --> exiting"; exit $rc; fi
novoalign -o SAM "@RG\tID:${READGROUP}\tSM:${sample}\tLB:${lib}" -d "$INDEXNAME" -f "$readu" -F "$FQFMT" -R "$MAPDIFF" -r "$REPEAT" -i PE "$AVGSZ","$SDSZ" -x "$GAPEXTEND" -l "$MINBASE" > "${OUTDIR}/NovoUnpaired.sam"
rc=$?; if [[ $rc != 0 ]]; then echo "Call to Novoalign failed on unpaired library $fullname with return code $rc --> exiting"; exit $rc; fi
else
	if [ "$MAXSCORE" -le 0 ]
	then
	>&2 echo "Error: Maximum alignment score must be positive --> exiting"
	exit 1
	else
	>&2 echo "Maximum alignment score: $MAXSCORE"
	fi

novoalign -o SAM "@RG\tID:${READGROUP}\tSM:${sample}\tLB:${lib}" -d "$INDEXNAME" -f "$read1" "$read2" -F "$FQFMT" -t "$MAXSCORE" -R "$MAPDIFF" -r "$REPEAT" -i PE "$AVGSZ","$SDSZ" -x "$GAPEXTEND" -l "$MINBASE" > "${OUTDIR}/NovoPaired.sam"
rc=$?; if [[ $rc != 0 ]]; then echo "Call to Novoalign failed on paired library $fullname with return code $rc --> exiting"; exit $rc; fi
novoalign -o SAM "@RG\tID:${READGROUP}\tSM:${sample}\tLB:${lib}" -d "$INDEXNAME" -f "$readu" -F "$FQFMT" -t "$MAXSCORE" -R "$MAPDIFF" -r "$REPEAT" -i PE "$AVGSZ","$SDSZ" -x "$GAPEXTEND" -l "$MINBASE" > "${OUTDIR}/NovoUnpaired.sam"
rc=$?; if [[ $rc != 0 ]]; then echo "Call to Novoalign failed on unpaired library $fullname with return code $rc --> exiting"; exit $rc; fi
fi

# extract only reads that aligned to reference
# ZS:Z:NM = No alignment was found
# ZS:Z:QC = Read failed quality checks and was not aligned

grep -v 'ZS:Z:NM\|ZS:Z:QC' "${OUTDIR}/NovoPaired.sam" > "${OUTDIR}/NovoPaired_aligned_only.sam"
rc=$?; if [[ $rc != 0 ]]; then echo "Failed to grep only aligned reads for paired library $fullname --> exiting"; exit $rc; fi
grep -v 'ZS:Z:NM\|ZS:Z:QC' "${OUTDIR}/NovoUnpaired.sam" > "${OUTDIR}/NovoUnpaired_aligned_only.sam"
rc=$?; if [[ $rc != 0 ]]; then echo "Failed to grep only aligned reads for unpaired library $fullname --> exiting"; exit $rc; fi

# convert SAM to BAM

samtools view -b -q $MINMAPQ "${OUTDIR}/NovoPaired_aligned_only.sam" > "${OUTDIR}/NovoPaired_aligned_only.bam"
rc=$?; if [[ $rc != 0 ]]; then echo "Call to SAMtools 'view' failed with return code $rc for paired library $fullname --> exiting"; exit $rc; fi
samtools view -b -q $MINMAPQ "${OUTDIR}/NovoUnpaired_aligned_only.sam" > "${OUTDIR}/NovoUnpaired_aligned_only.bam"
rc=$?; if [[ $rc != 0 ]]; then echo "Call to SAMtools 'view' failed with return code $rc for unpaired library $fullname --> exiting"; exit $rc; fi

# combined paired and unpaired BAMs

samtools merge "${OUTDIR}/NovoCombined.bam" "${OUTDIR}/NovoPaired_aligned_only.bam" "${OUTDIR}/NovoUnpaired_aligned_only.bam"
rc=$?; if [[ $rc != 0 ]]; then echo "Call to SAMtools 'merge' failed with return code $rc for library $fullname --> exiting"; exit $rc; fi

# sort BAM file

outbam="${OUTDIR}/${fullname}_sorted.bam"
samtools sort -o $outbam "${OUTDIR}/NovoCombined.bam"
rc=$?; if [[ $rc != 0 ]]; then echo "Call to SAMtools 'sort' failed with return code $rc for library $fullname --> exiting"; exit $rc; fi

# index BAM file

samtools index "${outbam}"
rc=$?; if [[ $rc != 0 ]]; then echo "Call to SAMtools 'index' failed with return code $rc --> exiting"; exit $rc; fi

# delete intermediate files

rm "${OUTDIR}/NovoPaired.sam"
rm "${OUTDIR}/NovoUnpaired.sam"
rm "${OUTDIR}/NovoPaired_aligned_only."*am
rm "${OUTDIR}/NovoUnpaired_aligned_only."*am
rm "${OUTDIR}/NovoCombined.bam"

done

>&2 echo "Finished!"

exit
