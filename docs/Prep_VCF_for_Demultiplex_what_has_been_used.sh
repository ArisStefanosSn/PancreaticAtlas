#!/bin/bash

#SBATCH -N 1
#SBATCH -n 20
#SBATCH --time 12:00:00
#SBATCH --mem=160GB
#SBATCH --mail-user=redwan.bhuiyan@jax.org
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute -q batch
#SBATCH --job-name=Yale.moving.to.hg38

cd $SLURM_SUBMIT_DIR

module load singularity
module load gcc

export PATH=$PATH:/path/to/Programs/plink-1.9
export PATH=$PATH:/path/to/Programs/bcftools-1.11

#Make binary PLINK files
plink --make-bed --file Stitzel_GS_010721 --out Stitzel_GS_010721

#Make freq files
plink --freq \
 --bfile Stitzel_GS_010721 \
 --out Stitzel_GS_010721

#Check to see if the bim and frq files have been properly created
perl /path/to/Programs/HRC-1000G-check-bim-NoReadKey.pl -b Stitzel_GS_010721.bim -f Stitzel_GS_010721.frq -r /path/to/Programs/PASS.Variantsbravo-dbsnp-all.tab -h

#Making the appropriate adjustments to the binary files (Run as is)
sh Run-plink.sh

#Combining all chromosomal vcf files together
bcftools concat Stitzel_GS_010721-updated-chr*.vcf -o Stitzel_GS_010721-updated-combined.vcf

#Sort VCF files
bcftools sort -m 30G Stitzel_GS_010721-updated-combined.vcf -o Stitzel_GS_010721-updated-combined-sorted.vcf

#Add chr prefix in front of chromosome number
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' Stitzel_GS_010721-updated-combined-sorted.vcf > Stitzel_GS_010721-updated-combined-sorted-with-chr.vcf

#Liftover from hg19 to hg38
java -jar /path/to/Programs/picard.jar LiftoverVcf \
 I=Stitzel_GS_010721-updated-combined-sorted-with-chr.vcf \
 O=Stitzel_GS_010721-updated-combined-sorted-hg38.vcf \
 CHAIN=/path/to/Reference/Genomes/hg19ToHg38.over.chain \
 REJECT=rejected_variants.vcf \
 R=/path/to/Reference/Genomes/hg38/hg38.fa

#Filter to get the proper samples
bcftools view --threads 20 -S Samples_filter.txt -o Stitzel_GS_010721-updated-combined-sorted-hg38-filtered.vcf Stitzel_GS_010721-updated-combined-sorted-hg38.vcf

#Rename the VCF file sample names to match Islet Numbers
bcftools reheader -s Samples_reheader.txt -o Siddhi_and_Makis.vcf Stitzel_GS_010721-updated-combined-sorted-hg38-filtered.vcf
