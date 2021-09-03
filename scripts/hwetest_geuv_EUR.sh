#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l mem=24gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -t 1-21
#PBS -N HWE 
#PBS -j oe
#PBS -o /home/vitor/heritability-hla-expression/log/$PBS_JOBID

CHR=$PBS_ARRAYID
SAMPLES=/home/vitor/heritability-hla-expression/data/ids_geuvadis_EUR.txt
VCFIN=/raid/genevol/vcf_1000G/phase3_20130502_grch38positions/ALL.chr${CHR}_GRCh38.genotypes.20170504.vcf.gz
OUT=/scratch/vitor/chr${CHR}_geuv_eur

bcftools view --samples-file $SAMPLES --force-samples $VCFIN |\
    bcftools view --genotype ^miss - |\
    bcftools view --min-af 0.01:minor - |\
    bcftools norm -m +both - |\
    vcftools --vcf - --hardy --out $OUT

