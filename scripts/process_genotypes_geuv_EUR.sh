#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l mem=24gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -t 1-22
#PBS -N processVCF
#PBS -j oe
#PBS -o /home/vitor/heritability-hla-expression/log/$PBS_JOBIB

CHR=$PBS_ARRAYID
SAMPLES=/home/vitor/heritability-hla-expression/data/ids_geuvadis_EUR.txt
VCFIN=/raid/genevol/vcf_1000G/phase3_20130502_grch38positions/ALL.chr${CHR}_GRCh38.genotypes.20170504.vcf.gz
VCFOUT=/scratch/vitor/chr${CHR}_geuv_eur.vcf
GDS=/scratch/vitor/chr${CHR}_geuv_eur.gds 
SCRIPTS=/home/vitor/heritability-hla-expression/scripts

bcftools view --samples-file $SAMPLES --force-samples $VCFIN |\
    bcftools view --genotype ^miss - |\
    bcftools view --min-af 0.01:minor - |\
    vcftools --vcf - --hwe 0.00000005 --recode --recode-INFO-all --stdout |\
    bcftools norm -m +both -o $VCFOUT -

Rscript $SCRIPTS/vcf2gds.R $VCFOUT $GDS 
rm $VCFOUT

Rscript $SCRIPTS/prune.R $GDS

