#!/bin/bash

# This runs kallisto quant for all fastq.gz files (must be run from directory with series).

SPECIES=$1
PAIRED=$2

sudo chmod -R a+rwX kallisto/
sudo chmod -R a+rwX logs/

if [[ $PAIRED = "" || $SPECIES = "" ]]; then
    echo "usage kallisto_dirrun.sh vXX p/s [optional] xxx/yyy/zzz_"
elif [[ $PAIRED = "p" ]]; then
    for FQ in fastqs/*.R1.fastq.gz; do
        while [ $(jobs | wc -l) -ge 4 ] ; do sleep 1 ; done
        ~/various/kallisto_paired.sh ${FQ%%.R1.fastq.gz} $SPECIES &
    done
elif [[ $PAIRED = "s" ]]; then
    for FQ in fastqs/*.fastq.gz; do
        while [ $(jobs | wc -l) -ge 4 ] ; do sleep 1 ; done
        ~/various/kallisto_single.sh ${FQ%%.fastq.gz} $SPECIES &
    done    
fi

for TSV in kallisto/raw/*_kallisto_out; do 
    TRUE=${TSV##kallisto/raw/}
    cp $TSV/abundance.tsv kallisto/tsvs/${TRUE%%_kallisto_out}.tsv; done
    
echo All over here

    
