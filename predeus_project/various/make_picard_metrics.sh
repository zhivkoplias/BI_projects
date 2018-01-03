#!/bin/bash

for i in bams/*.bam 
do 
    i=${i#bams/}
    TAG=${i%%.bam}
    sudo /stats/star_stat_hm38G.sh bams STAR_logs $TAG > stats/${TAG}.picard.log &
done
wait

echo -e "Sample\tN_reads\tPct_mapped\tPct_mapped_1loc\tPct_unmapped\tPct_rRNA\tPct_coding\tPct_UTR\tPct_intronic\tPct_intergenic\tJunctions\tInsertion_rate\tDeletion_rate\tPct_NC_junctions\tDel_av_length\tIns_av_length" > folder.stats
cat *rnastat >> folder.stats


