#!/bin/bash

#SBATCH -c 4
#SBATCH -N 1
#SBATCH -t 0-02:00
#SBATCH -p priority
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

cd /n/groups/cbdm-db/dam41/FC_06785

computeMatrix reference-point \
	-S  alignments/bw/MEC_Grhl1.merged.CPM.bigwig \
	-R 	deeptools/bed/shared_idr_peaks_minusIgG_20220215 \
        deeptools/bed/scatac_peaks/unique_6_idr_peaks \
        deeptools/bed/scatac_peaks/unique_1_idr_peaks \
		deeptools/bed/scatac_peaks/unique_2_idr_peaks \
		deeptools/bed/scatac_peaks/unique_3_idr_peaks \
		deeptools/bed/scatac_peaks/unique_4_idr_peaks \
		deeptools/bed/scatac_peaks/unique_5_idr_peaks \
		deeptools/bed/scatac_peaks/unique_10_idr_peaks \
	-out deeptools/matrices/Grhl1_merged_scatac_peaks_5.matrix \
	--samplesLabel Grhl1 \
	--referencePoint center \
	-b 5000 \
	-a 5000	\
	-p 8 \
	-bl /n/groups/cbdm-db/dam41/FC_06785/reference/mm10.blacklist.bed \
	--verbose
	
wait	

plotProfile -m deeptools/matrices/Grhl1_merged_scatac_peaks_5.matrix \
	-o deeptools/pdf/Grhl1_merged_scatac_peaks_profile_4.pdf \
	--verbose --perGroup --regionsLabel Shared C6 C1 C2 C3 C4 C5 C10 \
	--yMin 0 --yMax 9
plotHeatmap -m deeptools/matrices/Grhl1_merged_scatac_peaks_5.matrix \
	-o deeptools/pdf/Grhl1_merged_scatac_peaks_heatmap_4.pdf \
	--verbose --colorMap Blues --perGroup \
	--missingDataColor 1 --regionsLabel Shared C6 C1 C2 C3 C4 C5 C10 \
	--yMin 0.1 --yMax 9 --zMin 0.1 --zMax 8

computeMatrix reference-point \
	-S  alignments/bw/MEC_Hnf4a.merged.CPM.bigwig \
	-R 	deeptools/bed/shared_idr_peaks_minusIgG_20220215 \
        deeptools/bed/scatac_peaks/unique_6_idr_peaks \
        deeptools/bed/scatac_peaks/unique_1_idr_peaks \
		deeptools/bed/scatac_peaks/unique_2_idr_peaks \
		deeptools/bed/scatac_peaks/unique_3_idr_peaks \
		deeptools/bed/scatac_peaks/unique_4_idr_peaks \
		deeptools/bed/scatac_peaks/unique_5_idr_peaks \
		deeptools/bed/scatac_peaks/unique_10_idr_peaks \
	-out deeptools/matrices/Hnf4a_merged_scatac_peaks_5.matrix \
	--samplesLabel Hnf4a \
	--referencePoint center \
	-b 5000 \
	-a 5000	\
	-p 8 \
	-bl /n/groups/cbdm-db/dam41/FC_06785/reference/mm10.blacklist.bed \
	--verbose
	
wait	

plotProfile -m deeptools/matrices/Hnf4a_merged_scatac_peaks_5.matrix \
	-o deeptools/pdf/Hnf4a_merged_scatac_peaks_profile_4.pdf --verbose \
	--perGroup --regionsLabel Shared C6 C1 C2 C3 C4 C5 C10 \
	--yMin 0.1 --yMax 0.9
plotHeatmap -m deeptools/matrices/Hnf4a_merged_scatac_peaks_5.matrix \
	-o deeptools/pdf/Hnf4a_merged_scatac_peaks_heatmap_4.pdf --verbose \
	--colorMap Reds --perGroup --missingDataColor 1 --regionsLabel Shared C6 C1 C2 C3 C4 C5 C10 \
	--yMin 0.1 --yMax 0.9 --zMin 0.1 --zMax 0.8

computeMatrix reference-point \
	-S  alignments/bw/MEC_IgG.merged.CPM.bigwig \
	-R 	deeptools/bed/shared_idr_peaks_minusIgG_20220215 \
        deeptools/bed/scatac_peaks/unique_6_idr_peaks \
        deeptools/bed/scatac_peaks/unique_1_idr_peaks \
		deeptools/bed/scatac_peaks/unique_2_idr_peaks \
		deeptools/bed/scatac_peaks/unique_3_idr_peaks \
		deeptools/bed/scatac_peaks/unique_4_idr_peaks \
		deeptools/bed/scatac_peaks/unique_5_idr_peaks \
		deeptools/bed/scatac_peaks/unique_10_idr_peaks \
	-out deeptools/matrices/IgG_merged_scatac_peaks_5.matrix \
	--samplesLabel IgG \
	--referencePoint center \
	-b 5000 \
	-a 5000	\
	-p 8 \
	-bl /n/groups/cbdm-db/dam41/FC_06785/reference/mm10.blacklist.bed \
	--verbose
	
wait	

plotProfile -m deeptools/matrices/IgG_merged_scatac_peaks_5.matrix \
	-o deeptools/pdf/IgG_merged_scatac_peaks_profile_4.pdf --verbose \
	--perGroup --regionsLabel Shared C6 C1 C2 C3 C4 C5 C10 \
	--yMin 0.1 --yMax 10
plotHeatmap -m deeptools/matrices/IgG_merged_scatac_peaks_5.matrix \
	-o deeptools/pdf/IgG_merged_scatac_peaks_heatmap_4.pdf --verbose \
	--colorMap Greys --perGroup --missingDataColor 1 --regionsLabel Shared C6 C1 C2 C3 C4 C5 C10 \
	--yMin 0.1 --yMax 50 --zMin 0.1 --zMax 10

computeMatrix reference-point \
	-S  /n/groups/cbdm-db/dam41/FC_PostAireTF/alignments/bw/MEC_H3K27ac.merged.CPM.bigwig \
	-R 	deeptools/bed/shared_idr_peaks_minusIgG_20220215 \
        deeptools/bed/scatac_peaks/unique_6_idr_peaks \
        deeptools/bed/scatac_peaks/unique_1_idr_peaks \
		deeptools/bed/scatac_peaks/unique_2_idr_peaks \
		deeptools/bed/scatac_peaks/unique_3_idr_peaks \
		deeptools/bed/scatac_peaks/unique_4_idr_peaks \
		deeptools/bed/scatac_peaks/unique_5_idr_peaks \
		deeptools/bed/scatac_peaks/unique_10_idr_peaks \
	-out deeptools/matrices/H3K27ac_merged_scatac_peaks_5.matrix \
	--samplesLabel H3K27ac \
	--referencePoint center \
	-b 5000 \
	-a 5000 \
	-p 8 \
	-bl /n/groups/cbdm-db/dam41/FC_06785/reference/mm10.blacklist.bed \
	--verbose
	
wait	

plotProfile -m deeptools/matrices/H3K27ac_merged_scatac_peaks_5.matrix \
	-o deeptools/pdf/H3K27ac_merged_scatac_peaks_profile_4.pdf --verbose \
	--perGroup --regionsLabel Shared C6 C1 C2 C3 C4 C5 C10 \
	--yMin 0.1 --yMax 8
plotHeatmap -m deeptools/matrices/H3K27ac_merged_scatac_peaks_5.matrix \
	-o deeptools/pdf/H3K27ac_merged_scatac_peaks_heatmap_4.pdf --verbose \
	--colorMap Purples --perGroup --missingDataColor 1 --regionsLabel Shared C6 C1 C2 C3 C4 C5 C10 \
	--yMin 0.1 --yMax 8 --zMin 0.1 --zMax 2

computeMatrix reference-point \
	-S  alignments/bw/MEC_Aire.merged.CPM.bigwig \
	-R 	deeptools/bed/shared_idr_peaks_minusIgG_20220215 \
        deeptools/bed/scatac_peaks/unique_6_idr_peaks \
        deeptools/bed/scatac_peaks/unique_1_idr_peaks \
		deeptools/bed/scatac_peaks/unique_2_idr_peaks \
		deeptools/bed/scatac_peaks/unique_3_idr_peaks \
		deeptools/bed/scatac_peaks/unique_4_idr_peaks \
		deeptools/bed/scatac_peaks/unique_5_idr_peaks \
		deeptools/bed/scatac_peaks/unique_10_idr_peaks \
	-out deeptools/matrices/Aire_merged_scatac_peaks_5.matrix \
	--samplesLabel Aire \
	--referencePoint center \
	-b 5000 \
	-a 5000 \
	-p 8 \
	-bl /n/groups/cbdm-db/dam41/FC_06785/reference/mm10.blacklist.bed \
	--verbose
	
wait	

plotProfile -m deeptools/matrices/Aire_merged_scatac_peaks_5.matrix \
	-o deeptools/pdf/Aire_merged_scatac_peaks_profile_4.pdf --verbose \
	--perGroup --regionsLabel Shared C6 C1 C2 C3 C4 C5 C10 \
	--yMin 0.01 --yMax 0.6
plotHeatmap -m deeptools/matrices/Aire_merged_scatac_peaks_5.matrix \
	-o deeptools/pdf/Aire_merged_scatac_peaks_heatmap_4.pdf --verbose \
	--colorMap Greens --perGroup --missingDataColor 1 --regionsLabel Shared C6 C1 C2 C3 C4 C5 C10 \
	--yMin 0.01 --yMax 0.6 --zMin 0.01 --zMax 0.3

cd /n/groups/cbdm-db/dam41/FC_06921	

computeMatrix reference-point \
	-S  alignments/bw/MEC_Pou2f3_all.merged.CPM.bigwig \
	-R 	deeptools/bed/shared_idr_peaks_minusIgG_20220215 \
		deeptools/bed/scatac_peaks/unique_6_idr_peaks \
        deeptools/bed/scatac_peaks/unique_1_idr_peaks \
		deeptools/bed/scatac_peaks/unique_2_idr_peaks \
		deeptools/bed/scatac_peaks/unique_3_idr_peaks \
		deeptools/bed/scatac_peaks/unique_4_idr_peaks \
		deeptools/bed/scatac_peaks/unique_5_idr_peaks \
		deeptools/bed/scatac_peaks/unique_10_idr_peaks \
	-out deeptools/matrices/Pou2f3_all_merged_scatac_peaks3.matrix \
	--samplesLabel Pou2f3_all \
	--referencePoint center \
	-b 5000 \
	-a 5000	\
	-p 8 \
	--verbose
	
wait	

plotProfile -m deeptools/matrices/Pou2f3_all_merged_scatac_peaks3.matrix \
	-o deeptools/pdf/Pou2f3_all_merged_scatac_peaks_profile_4.pdf \
	--verbose --perGroup --regionsLabel Shared C6 C1 C2 C3 C4 C5 C10 \
	--yMin 0 --yMax 2.5
plotHeatmap -m deeptools/matrices/Pou2f3_all_merged_scatac_peaks3.matrix \
	-o deeptools/pdf/Pou2f3_all_merged_scatac_peaks_heatmap_4.pdf \
	--verbose --colorMap Oranges --perGroup \
	--missingDataColor 1 --regionsLabel Shared C6 C1 C2 C3 C4 C5 C10 \
	--yMin 0 --yMax 2.5 --zMin 0 --zMax 2.5
