#!/bin/bash
#PBS -l nodes=1:ppn=1,mem=20gb
#PBS -l walltime=4:00:00

module load bedtools/2.25.0

cd $PBS_O_WORKDIR

chr_size=/nfs/rprdata/Anthony/data/hg38_annotations/hg38_ref_chrom_order.txt
bedtranscript=/nfs/rprdata/Anthony/data/hg38_annotations/hg38_transcriptome.anno.no_head.bam_sort.tsv.gz
dataFolder=../Bam

bedtools coverage -sorted -counts -split -s -sorted -a ${bedtranscript} \
            -b ${dataFolder}/${Sample}_clean.bam -g ${chr_size} | gzip > \
             GeneCounts/${Sample}_counts.bed.gz

echo ${Sample} >> Finished.txt
