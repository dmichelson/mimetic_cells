#!/bin/bash

#SBATCH -c 4
#SBATCH -N 1
#SBATCH -t 0-02:00
#SBATCH -p short
#SBATCH --mem=4G
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
module load samtools/1.3.1
#module load bedtools/2.26.0
#module load homer/4.9
#module load multiqc/1.5
#module load trimmomatic/0.36
module load R/3.4.1
module load deeptools/3.0.2
#module load star/2.5.4a
#module load rsem/1.3.0
#module load picard/2.8.0
module load macs2/2.1.1.20160309

Sample=$1
Genome=$2
WorkingDir=$3

cache="/n/scratch3/users/d/dam41"

cd $WorkingDir
pwd

echo "Sample#=" $Sample "; Genome=" $Genome "; Working directory=" $WorkingDir

if [[ "$Sample" == "MEC_IgG_CT" ]]; then
	echo "Running samtools merge IgG ... "
	samtools merge -@ 4 $cache/${Sample}.merged.bam alignments/${Sample}_R1.uniq.bam alignments/${Sample}_R2.uniq.bam
fi

if [[ "$Sample" == "MEC_H3K27ac_CT" ]]; then
	echo "Running samtools merge H3K27ac ... "
	samtools merge -@ 4 $cache/${Sample}.merged.bam alignments/${Sample}_R1.uniq.bam alignments/${Sample}_R2.uniq.bam
fi

if [[ "$Sample" == "MEC_Grhl1" ]]; then
	echo "Running samtools merge Grhl1 ... "
	samtools merge -@ 4 $cache/${Sample}.merged.bam alignments/${Sample}_R2.uniq.bam alignments/${Sample}_R3.uniq.bam alignments/${Sample}_R4.uniq.bam
fi

if [[ "$Sample" == "MEC_Hnf4a" ]]; then
	echo "Running samtools merge Hnf4a ... "
	samtools merge -@ 4 $cache/${Sample}.merged.bam alignments/${Sample}_R1.uniq.bam alignments/${Sample}_R2.uniq.bam alignments/${Sample}_R3.uniq.bam alignments/${Sample}_R4.uniq.bam
fi

if [[ "$Sample" == "MEC_Pou2f3_all" ]]; then
	echo "Running samtools merge Pou2f3 ... "
	samtools merge -@ 4 $cache/${Sample}.merged.bam \
	alignments/MEC_Pou2f3_fix_R1.uniq.bam \
	alignments/MEC_Pou2f3_fix_R2.uniq.bam \
	alignments/MEC_Pou2f3_fix_R3.uniq.bam \
	alignments/MEC_Pou2f3_fix_R4.uniq.bam \
	alignments/MEC_Pou2f3_nofix_R1.uniq.bam \
	alignments/MEC_Pou2f3_nofix_R2.uniq.bam \
	alignments/MEC_Pou2f3_nofix_R3.uniq.bam \
	alignments/MEC_Pou2f3_nofix_R4.uniq.bam
fi

echo "Running samtools sort ... "
samtools view -@ 4 -bS $cache/${Sample}.merged.bam | samtools sort -@ 4 > alignments/bam_merge/${Sample}.merged.bam

echo "Running samtools index ... "
samtools index alignments/bam_merge/${Sample}.merged.bam

# echo "Running samtools ecoli merge ... "
# samtools merge -f -@ 4 $cache/${Sample}.ecoli.merged.bt2out.sam $cache/${Sample}_R1.ecoli.bt2out.sam $cache/${Sample}_R2.ecoli.bt2out.sam $cache/${Sample}_R3.ecoli.bt2out.sam $cache/${Sample}_R4.ecoli.bt2out.sam

echo "Running samtools ecoli merge ... "
samtools merge -f -@ 4 $cache/${Sample}.ecoli.merged.bt2out.sam \
	$cache/MEC_Pou2f3_fix_R1.ecoli.bt2out.sam \
	$cache/MEC_Pou2f3_fix_R2.ecoli.bt2out.sam \
	$cache/MEC_Pou2f3_fix_R3.ecoli.bt2out.sam \
	$cache/MEC_Pou2f3_fix_R4.ecoli.bt2out.sam \
	$cache/MEC_Pou2f3_nofix_R1.ecoli.bt2out.sam \
	$cache/MEC_Pou2f3_nofix_R2.ecoli.bt2out.sam \
	$cache/MEC_Pou2f3_nofix_R3.ecoli.bt2out.sam \
	$cache/MEC_Pou2f3_nofix_R4.ecoli.bt2out.sam

echo "Running samtools ecoli sort ... "
samtools view -@ 4 -bS $cache/${Sample}.ecoli.merged.bt2out.sam | samtools sort -@ 4 > $cache/${Sample}.ecoli.merged.sorted.sam

echo "Running samtools ecoli index ... "
samtools index $cache/${Sample}.ecoli.merged.sorted.sam

seqDepthDouble=`samtools view -F 0x04 $cache/${Sample}.ecoli.merged.sorted.sam | wc -l`
seqDepth=$((seqDepthDouble/2))
echo $seqDepth > alignments/qc/${Sample}.ecoli.merged.seqDepth

seqDepth=$(cat alignments/qc/${Sample}.ecoli.merged.seqDepth)
scaleFactor=`echo "1000000 / $seqDepth" | bc -l`
echo "seqDepth for " $Sample " is " $seqDepth "; scalingFactor is " $scaleFactor

echo "Making bedgraph file ..."
bamCoverage \
	-b alignments/bam_merge/${Sample}.merged.bam \
	-o alignments/bdg/${Sample}.merged.scaleFactor.bedgraph \
	-of bedgraph \
	-e \
	-p 4 \
	--scaleFactor $scaleFactor

echo "Making bigwig file ..."
bamCoverage \
	-b alignments/bam_merge/${Sample}.merged.bam \
	-o alignments/bw/${Sample}.merged.scaleFactor.bigwig \
	-of bigwig \
	-e \
	-p 4 \
	--scaleFactor $scaleFactor

echo "Making bedgraph file ..."
bamCoverage \
	-b alignments/bam_merge/${Sample}.merged.bam \
	-o alignments/bdg/${Sample}.merged.CPM.bedgraph \
	-of bedgraph \
	-p 4 \
	--normalizeUsing CPM \
    # -e

echo "Making bigwig file ..."
bamCoverage \
	-b alignments/bam_merge/${Sample}.merged.bam \
	-o alignments/bw/${Sample}.merged.CPM.bigwig \
	-of bigwig \
	-p 4 \
	--normalizeUsing CPM \
    # -e 

echo "Running macs2 ..."
macs2 callpeak -t alignments/bam_merge/${Sample}.merged.bam -f BAMPE -g mm -q 0.05 -n $Sample --outdir peaks/macs/merge/bampe
macs2 callpeak -t alignments/bam_merge/${Sample}.merged.bam -f BAM -g mm -q 0.05 --shift 100 --extsize 200 -n $Sample --outdir peaks/macs/merge/shift

### call peaks with seacr
seacr="SEACR/SEACR_1.3.sh"
# controlSample=MEC_IgG

# bash $seacr \
# 	alignments/bdg/${Sample}.merged.scaleFactor.bedgraph \
#     alignments/bdg/${controlSample}.merged.scaleFactor.bedgraph \
#     non stringent \
# 	peaks/seacr/merge/${Sample}_merged_seacr_control.peaks

bash $seacr \
	alignments/bdg/${Sample}.merged.scaleFactor.bedgraph \
	0.01 non stringent \
	peaks/seacr/merge/${Sample}_merged_seacr_top0.01.peaks

echo "Mapping stats... we found this many merged reads... "
samtools view alignments/bam_merge/${Sample}.merged.bam | wc -l
