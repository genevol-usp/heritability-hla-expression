#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l mem=24gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -t 1-22
#PBS -N processVCF
#PBS -j oe
#PBS -o /home/vitor/heritability-hla-expression/log/$PBS_JOBID

CHR=$PBS_ARRAYID
VCFIN=/raid/genevol/vcf_1000G/phase3_20130502_grch38positions/ALL.chr${CHR}_GRCh38.genotypes.20170504.vcf.gz
VCFOUT=/scratch/vitor/chr${CHR}.vcf
GDS=/scratch/vitor/chr${CHR}.gds 
SCRIPTS=/home/vitor/heritability-hla-expression/scripts

bcftools view --genotype ^miss $VCFIN |\
    bcftools view --min-af 0.01:minor - |\
    vcftools --vcf - --hwe 0.00000005 --recode --recode-INFO-all --stdout |\
    bcftools norm -m +both -o $VCFOUT -

Rscript $SCRIPTS/vcf2gds.R $VCFOUT $GDS 
rm $VCFOUT

Rscript $SCRIPTS/prune.R $GDS

