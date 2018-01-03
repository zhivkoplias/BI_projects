#!/bin/bash

# This downloads whole series using tab-delimited format file

SAMPLES=$1

links=`awk -F "\t" '{ if (NR > 1)  print $31 }' $SAMPLES`

for link in $links
do
    wget $link
done

awk -F "\t" '{ if (NR > 1)  print $30 }' $SAMPLES > errs
awk -F "\t" '{ if (NR > 1)  print $1 }' $SAMPLES > samples
paste errs samples | sort -u > SAMPLES
rm errs samples
