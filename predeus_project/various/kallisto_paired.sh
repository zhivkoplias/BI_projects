#!/bin/bash 

# Running kallisto for pe reads

PTH=$1
SPECIES=$2

TAG=${PTH##fastqs/}
echo $TAG
SAMPLES="$(cat SAMPLES)"
if [[ $SAMPLES == *$TAG* ]]; then
	gunzip -c ${PTH}.R1.fastq.gz > ${TAG}.R1.fastq
	gunzip -c ${PTH}.R2.fastq.gz > ${TAG}.R2.fastq
	kallisto quant -i /media/DISK1/reference/Kallisto_indexes/gencode.${SPECIES}.index -o kallisto/raw/${TAG}_kallisto_out -t 4 ${TAG}.R1.fastq ${TAG}.R2.fastq &> logs/${TAG}.kallisto.log
	rm ${TAG}.R1.fastq
	rm ${TAG}.R2.fastq
else
	echo Skipping...
fi

