#!/bin/bash

#SBATCH -c 4
#SBATCH -N 1
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem-per-cpu=2G
#SBATCH -o logs/dam41_%j.out
#SBATCH -e logs/dam41_%j.err

module load gcc/6.2.0
module load python/2.7.12
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
module load deeptools/3.0.2
#module load star/2.5.4a
#module load rsem/1.3.0
#module load picard/2.8.0

cd /n/groups/cbdm-db/dam41/FC_06921

Sample=$1
Bigwig=$2
Region1=$3
Region2=$4

computeMatrix reference-point \
	-S $Bigwig \
	-R $Region1 \
	   $Region2 \
	-out deeptools/matrices/${Sample}_plusMinus_aire.matrix \
	--samplesLabel $Sample \
	--referencePoint center \
	-b 2500 \
	-a 2500	\
	-bl reference/mm10.blacklist.bed \
	-p 4 \
	--verbose
	
wait	

plotProfile -m deeptools/matrices/${Sample}_plusMinus_aire.matrix -o deeptools/pdf/plusMinusAire/${Sample}_plusMinus_aire_prof.pdf \
	--verbose --perGroup --regionsLabel MinusAire PlusAire \
	--yMin 0 --yMax 4
plotHeatmap -m deeptools/matrices/${Sample}_plusMinus_aire.matrix -o deeptools/pdf/plusMinusAire/${Sample}_plusMinus_aire_heat.pdf \
	--verbose --colorMap Reds --perGroup --missingDataColor 1 --regionsLabel MinusAire PlusAire \
	--yMin 0 --yMax 4 --zMin 0 --zMax 4

