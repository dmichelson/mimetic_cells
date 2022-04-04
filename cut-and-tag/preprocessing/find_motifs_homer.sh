#!/bin/bash

#SBATCH -c 4
#SBATCH -N 1
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem-per-cpu=500M
#SBATCH -o logs/dam41_%j.out
#SBATCH -e logs/dam41_%j.err

module load gcc/6.2.0
# module load python/2.7.12
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
# module load macs2/2.1.1.20160309
export PATH=/n/scratch3/users/d/dam41/homer/bin:$PATH

Sample=$1
Genome=$2
WorkingDir=$3

cd $WorkingDir

mkdir peaks/homer/MEC_${Sample}_top1
mkdir peaks/homer/MEC_${Sample}_macs_bampe
mkdir peaks/homer/MEC_${Sample}_macs_shift

findMotifsGenome.pl peaks/seacr/merge/MEC_${Sample}_merged_seacr_top0.01.peaks.stringent.bed mm10r peaks/homer/MEC_${Sample}_top1 -p 4
findMotifsGenome.pl peaks/macs/merge/shift/MEC_${Sample}_peaks.narrowPeak mm10r peaks/homer/MEC_${Sample}_macs_shift -p 4
findMotifsGenome.pl peaks/macs/merge/bampe/MEC_${Sample}_peaks.narrowPeak mm10r peaks/homer/MEC_${Sample}_macs_bampe -p 4
