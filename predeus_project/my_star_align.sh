#!/bin/bash 

## can be used for both single and paired-end
## archived FASTQ is assumed

TAG=$1
RL=101
R1=$2
R2=$3
READS=""
WDIR=`pwd`

if [[ $TAG == "" || $R1 == "" ]]
then 
  echo "Please provide 1) output name (tag) 2) R1 (or R1/R2 for paired-end) read file name"
  exit 1
fi 

if [[ $R2 == "" ]]
then 
  echo "Processing alignment as single-end, using STAR index /media/DISK1/celllines/STAR/gencode_v23_101bp_101bp."
  READS=$WDIR/$R1
else 
  echo "Processing alignment as paired-end, using STAR index /media/DISK1/celllines/STAR/gencode_v23_101bp_101bp"
  READS="$WDIR/$R1 $WDIR/$R2"
fi

GENDIR="/media/DISK1/celllines/STAR/gencode_v23_101bp_101bp"


mkdir ${TAG}_STAR
cd ${TAG}_STAR
STAR --genomeDir $GENDIR --readFilesIn $READS --runThreadN 4 --readFilesCommand zcat --outFilterMultimapNmax 15 --outFilterMismatchNmax 6  --outSAMstrandField All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM 

mv Aligned.sortedByCoord.out.bam $TAG.bam
mv Aligned.toTranscriptome.out.bam $TAG.tr.bam 
mv Log.out $TAG.star_run.log 
mv Log.final.out $TAG.star_final.log

