#!/bin/bash

#SBATCH -c 16
#SBATCH -N 1
#SBATCH -t 1-12:00
#SBATCH -p medium
#SBATCH --mem-per-cpu=1G
#SBATCH -o logs/dam41_%j.out
#SBATCH -e logs/dam41_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=daniel_michelson@hms.harvard.edu

#module load gcc/6.2.0
#module load python/2.7.12
#module load java/jdk-1.8u112
#module load git/2.9.5
#module load bwa/0.7.8
#module load cuda/9.0
#module load sratoolkit/2.8.1
#module load fastqc/0.11.5
#module load bowtie2/2.2.9
#module load samtools/1.3.1
#module load bedtools/2.26.0
#module load homer/4.9
#module load multiqc/1.5
#module load trimmomatic/0.36
#module load R/3.4.1
#module load deeptools/3.0.2
#module load star/2.5.4a
#module load rsem/1.3.0
#module load picard/2.8.0
module load cellranger-ATAC/1.1.0

cd /n/groups/cbdm-db/dam41/mec_scatac/cellranger

cellranger-atac count \
    --id=$1 \
    --reference=/n/shared_db/mm10/uk/cellranger-ATAC/1.1.0/1.1.0/refdata-cellranger-atac-mm10-1.1.0 \
    --fastqs=/n/groups/cbdm-db/dam41/mec_scatac/data/200305_A00794_0197_BHHM5WDRXX/fastq/mTEC_scATAC_1/$1 \
    --force-cells 4000 \
    --localcores 16 \
    --localmem 16

# WT_Aire_MEChi1_10xATAC_5K_REA09522
# WT2_reprep_Aire_MEChi1_10xATAC_5K
# KO_Aire_MEChi1_10xATAC_5K_REA09522
# KO2_reprep_Aire_MEChi1_10xATAC_5K