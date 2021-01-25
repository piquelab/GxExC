#!/bin/bash
#PBS -l nodes=1:ppn=1,mem=100gb
#PBS -l walltime=4:00:00

module load bedtools/2.25.0

cd $PBS_O_WORKDIR


bedtranscript=/wsu/home/groups/piquelab/data/RefTranscriptome/Sorted_NoChr_ensGene.hg19.2014.bed.gz
dataFolder=/wsu/home/fx/fx78/fx7820/rprdata/Anthony/GxE_iPSC/Shallow/${Plate}/Bam

# mygenomefile.txt is a 2-column, tab separated file containing the chromosome names and number of bp
bedtools coverage -counts -sorted -split -s -a ${bedtranscript} \
            -b ${dataFolder}/${Sample}_clean.bam -g mygenomefile.txt | gzip > \
             GeneCounts/${Sample}_counts.bed.gz

echo ${Sample} >> Finished.txt
