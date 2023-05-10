#!/bin/bash

#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time 12:00:00
#SBATCH --mem=350GB
#SBATCH --mail-user=redwan.bhuiyan@jax.org
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute -q batch
#SBATCH --job-name=UVA.Pre.Imputation

cd $SLURM_SUBMIT_DIR

module load singularity
module load gcc

export PATH=$PATH:/projects/stitzel-lab/c-bhuiyr/Programs/
export PATH=$PATH:/projects/stitzel-lab/c-bhuiyr/Programs/bcftools-1.11/
export PATH=$PATH:/projects/stitzel-lab/c-bhuiyr/Programs/htslib-1.14/

#Set the number of chromosomes
for i in {1..22}
do
#Make bed files
plink --make-bed --bfile mega_islet_raw_nophenoA -chr ${i} --out BINARY.VCF/mega_islet_raw_nophenoA_chr${i}

#Make freq files
plink --freq --bfile BINARY.VCF/mega_islet_raw_nophenoA_chr${i} --out BINARY.VCF/mega_islet_raw_nophenoA_chr${i}

#Check to see if the bim and frq files have been properly created
perl /projects/stitzel-lab/c-bhuiyr/Programs/HRC-1000G-check-bim-NoReadKey.pl \ 
 -b BINARY.VCF/mega_islet_raw_nophenoA_chr${i}.bim \ 
 -f BINARY.VCF/mega_islet_raw_nophenoA_chr${i}.frq \ 
 -r /projects/stitzel-lab/c-bhuiyr/Imputation_Tools/PASS.Variantsbravo-dbsnp-all.tab -h 

#Making the appropriate adjustments to the binary files
sh BINARY.VCF/Run-plink.sh      # created by the perl, run as it is

#Sort VCF files
bcftools sort /projects/stitzel-lab/lawlon/Genotype/UVA_Islet_2019/plink/BINARY.VCF/mega_islet_raw_nophenoA_chr${i}-updated-chr${i}.vcf \ 
 -O z \ 
 -o /projects/stitzel-lab/lawlon/Genotype/UVA_Islet_2019/plink/BINARY.VCF/VCF/mega_islet_raw_nophenoA_chr${i}.vcf.gz 

#Check vcf files: make sure that the file is ready for use
python /projects/stitzel-lab/c-bhuiyr/Programs/checkVCF.py \ 
 -r /projects/stitzel-lab/c-bhuiyr/Genomes/human_g1k_v37.fasta \ 
 -o out \
 /projects/stitzel-lab/lawlon/Genotype/UVA_Islet_2019/plink/BINARY.VCF/VCF/mega_islet_raw_nophenoA_chr${i}.vcf.gz 

done

mkdir BINARY.VCF/WHOLE.BINARY/

#Getting a list of VCF files
ls /projects/stitzel-lab/lawlon/Genotype/UVA_Islet_2019/plink/BINARY.VCF/VCF/mega_islet_raw_nophenoA_chr*.vcf.gz > Chromosomal_File_List.txt

#Combining chromosomal vcf files together
bcftools concat -o BINARY.VCF/WHOLE.BINARY/mega_islet_raw_nophenoA-updated-combined.vcf \
 -f BINARY.VCF/VCF/Chromosomal_File_List.txt

#Sort VCF files
bcftools sort BINARY.VCF/WHOLE.BINARY/mega_islet_raw_nophenoA-updated-combined.vcf -o BINARY.VCF/WHOLE.BINARY/mega_islet_raw_nophenoA-updated-combined-sorted.vcf

#Liftover from hg19 to hg38
java -jar /projects/stitzel-lab/c-bhuiyr/Imputation_Tools/picard.jar LiftoverVcf \
 I=BINARY.VCF/WHOLE.BINARY/mega_islet_raw_nophenoA-updated-combined-sorted_with_chr.vcf \
 O=BINARY.VCF/WHOLE.BINARY/mega_islet_raw_nophenoA-updated-combined-sorted-hg38.vcf \
 CHAIN=/projects/stitzel-lab/c-bhuiyr/Imputation_Tools/hg19ToHg38.over.chain \
 REJECT=BINARY.VCF/WHOLE.BINARY/rejected_variants.vcf \
 R=/projects/stitzel-lab/c-bhuiyr/Genomes/hg38.fa

#Bgzipping and indexing the hg38-lifted VCF file
bgzip -@ 16 BINARY.VCF/WHOLE.BINARY/mega_islet_raw_nophenoA-updated-combined-sorted-hg38.vcf
tabix -p vcf BINARY.VCF/WHOLE.BINARY/mega_islet_raw_nophenoA-updated-combined-sorted-hg38.vcf.gz

#Subset the 4 islets needed: Islet61, Islet62, Islet63, and Islet68
bcftools view \
 --threads 16 \
 -s Islet69_Islet69,AS142N_AS142N,Islet29_Islet29,Islet50_Islet50 \
 -O z \
 -o BINARY.VCF/WHOLE.BINARY/mega_islet_raw_nophenoA-updated-combined-sorted-hg38-subsetted.vcf.gz \
 BINARY.VCF/WHOLE.BINARY/mega_islet_raw_nophenoA-updated-combined-sorted-hg38.vcf.gz

#Indexing the VCF file
tabix -p vcf BINARY.VCF/WHOLE.BINARY/mega_islet_raw_nophenoA-updated-combined-sorted-hg38-subsetted.vcf.gz

#Renaming the header to the actual islet samples
bcftools reheader \
 -s Sample_list.txt \
 --threads 16 \
 -o BINARY.VCF/WHOLE.BINARY/Final_VCF.vcf.gz \
 BINARY.VCF/WHOLE.BINARY/mega_islet_raw_nophenoA-updated-combined-sorted-hg38-subsetted.vcf.gz
