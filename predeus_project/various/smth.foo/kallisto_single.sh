#!/bin/bash

# Running kallisto for se reads

$PTH=$1
$SPECIES=$2

TAG=${PTH##fastqs/}
echo $TAG
gunzip -c fastqs/${TAG}.fastq.gz > ${TAG}.fastq
kallisto quant -i /media/DISK1/reference/Kallisto_indexes/gencode.${SPECIES}.index -o kallisto/raw/${TAG}_kallisto_out --single -l 190 -s 40 ${TAG}.fastq &> logs/${TAG}.kallisto.log
rm ${TAG}.fastq
