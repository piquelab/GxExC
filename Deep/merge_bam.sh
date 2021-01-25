#!/bin/bash
#PBS -l nodes=1:ppn=1,mem=50gb
#PBS -l walltime=24:00:00


module load hisat2/2.0.4
module load samtools/1.4

# Combine bam files and deduplicate. We only did this for the DCM plates since they were sequenced twice.

cd $PBS_O_WORKDIR

###Merge Reads###
samtools merge ${Sample}_merged.bam ../../P1_R2/Bam/${Sample}_clean.bam ../../P1_R3/Bam/${Sample}_clean.bam
samtools index ${Sample}_merged.bam
samtools view -c ${Sample}_merged.bam > ${Sample}_merged_count.txt


###Deduplication###
samtools rmdup ${Sample}_merged.bam ${Sample}_clean.bam
samtools index ${Sample}_clean.bam
samtools view -c ${Sample}_clean.bam > ${Sample}_clean_count.txt

echo ${Sample} >> finished.txt
