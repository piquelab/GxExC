#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=80gb
#PBS -l walltime=4:00:00

module load R

cd $PBS_O_WORKDIR

Rscript --vanilla QuASAR_prep_P1_R2.R ${Individual} ${Plate}

