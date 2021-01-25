#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=60gb
#PBS -l walltime=24:00:00


module load hisat2/2.0.4
module load samtools/1.4

cd $PBS_O_WORKDIR

# I do the following two commands because the FASTQ files have underscores instead of dashes like I like separating the plate from the barcode.
Plate=`echo $Sample | cut -d- -f1`
Barcode=`echo $Sample | cut -d- -f2`

#Update these if copied from another directory. This is for plate DLCL1R1
filePath=/wsu/home/groups/piquelab/OurData/Nextseq/IPSC_Deep/C202SC18120920/raw_data

genomeindex=/nfs/rprdata/Anthony/data/HISAT2_Index/grch38_snp_tran/genome_snp_tran


###Align Reads###
hisat2 -p 8 -x ${genomeindex} -1 ${filePath}/${Plate}_${Barcode}_1.fq.gz \
                              -2 ${filePath}/${Plate}_${Barcode}_2.fq.gz \
      2> ${Sample}_aligned.bam.e | samtools view -b1 - > ${Sample}_aligned.bam



###Sort Reads###
samtools sort -@ 8 -T tmp_${Sample}_aligned.bam -o ${Sample}_sorted.bam ${Sample}_aligned.bam
samtools index ${Sample}_sorted.bam
samtools view -c ${Sample}_sorted.bam > ${Sample}_sorted_count.txt


###Quality Filter###
samtools view -b1 -q10 ${Sample}_sorted.bam > ${Sample}_quality.bam
samtools index ${Sample}_quality.bam
samtools view -c ${Sample}_quality.bam > ${Sample}_quality_count.txt


###Deduplication###
samtools rmdup ${Sample}_quality.bam ${Sample}_clean.bam
samtools index ${Sample}_clean.bam
samtools view -c ${Sample}_clean.bam > ${Sample}_clean_count.txt

echo ${Sample} >> finished.txt
