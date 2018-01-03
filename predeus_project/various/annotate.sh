#!/bin/bash

ERRS=`awk '{print $2}' SAMPLES`
NAMES=`awk '{print $1}' SAMPLES`

tech=( $ERRS )
real=( $NAMES )

for i in `seq 0 ${#tech[@]}`
do
    mv ${tech[$i]}_1.fastq.gz ${real[$i]}.R1.fastq.gz
    mv ${tech[$i]}_2.fastq.gz ${real[$i]}.R2.fastq.gz
done
