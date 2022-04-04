#!/bin/bash

#SBATCH -c 4
#SBATCH -N 1
#SBATCH -t 0-00:30
#SBATCH -p short
#SBATCH --mem-per-cpu=1G
#SBATCH -o logs/dam41_%j.out
#SBATCH -e logs/dam41_%j.err

module load gcc/6.2.0
module load python/2.7.12
module load java/jdk-1.8u112
module load git/2.9.5
module load bwa/0.7.8
module load cuda/9.0
# module load sratoolkit/2.9.0
module load fastqc/0.11.5
module load fastx/0.0.13
#module load sickle/1.2
#module load cutadapt/1.14
module load bowtie2/2.2.9
module load samtools/1.3.1
module load bedtools/2.26.0
module load homer/4.9
module load multiqc/1.5
module load trimmomatic/0.36
module load R/3.5.1
module load deeptools/3.0.2
#module load star/2.5.4a
#module load rsem/1.3.0
module load picard/2.8.0
module load macs2/2.1.1.20160309

#input paired-end reads as: "./script.sh library sample genome /working/directory"

Library=$1
Sample=$2
Genome=$3
working_dir=$4
cache="/n/scratch3/users/d/dam41"

cd $working_dir
pwd

echo "Library=" $Library "; Sample#=" $Sample "; Genome=" $Genome "; Working directory=" $working_dir

bowtieDir="/n/groups/shared_databases/bowtie2_indexes/mm10"
genomeDir="/n/groups/shared_databases/bowtie2_indexes/mm10.fa"

echo "Merging fastq's"

bzcat data/${Library}_*_R1.fastq.bz2 > $cache/${Sample}_merged_R1.fastq
bzcat data/${Library}_*_R2.fastq.bz2 > $cache/${Sample}_merged_R2.fastq

echo "Doing fastqc... "
fastqc $cache/${Sample}_merged_R1.fastq --outdir data/qc
fastqc $cache/${Sample}_merged_R2.fastq --outdir data/qc

echo "Trimming low-quality reads by trimmomatic ..."
java -Xmx2048m -jar $TRIMMOMATIC/trimmomatic-0.36.jar PE \
	$cache/${Sample}_merged_R1.fastq \
	$cache/${Sample}_merged_R2.fastq \
	$cache/${Sample}_R1.trimmed.fastq \
	$cache/${Sample}_R1.unpaired.fastq \
	$cache/${Sample}_R2.trimmed.fastq \
	$cache/${Sample}_R2.unpaired.fastq \
	LEADING:3 \
	TRAILING:3 \
	SLIDINGWINDOW:4:15 \
	ILLUMINACLIP:reference/nextera_adapters.fasta:2:30:10

echo "Mapping reads to mm10..."
bowtie2 -q -p 4 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -x $bowtieDir \
	-1 $cache/${Sample}_R1.trimmed.fastq -2 $cache/${Sample}_R2.trimmed.fastq -S $cache/${Sample}.bt2out.sam \
	&> alignments/qc/${Sample}_bt2_summary.txt

echo "Mapping reads to Ecoli..."
eColiIndex="/n/groups/shared_databases/igenome/Escherichia_coli/NCBI/DH10B/Sequence/Bowtie2Index/genome"
bowtie2 -q -p 4 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 --no-overlap --no-dovetail -x $eColiIndex \
	-1 $cache/${Sample}_R1.trimmed.fastq -2 $cache/${Sample}_R2.trimmed.fastq -S $cache/${Sample}.ecoli.bt2out.sam \
	&> alignments/qc/${Sample}_ecoli_bt2_summary.txt

seqDepthDouble=`samtools view -F 0x04 $cache/${Sample}.ecoli.bt2out.sam | wc -l`
seqDepth=$((seqDepthDouble/2))
echo $seqDepth > alignments/qc/${Sample}.ecoli.seqDepth

echo "Removing multiple mapping reads... "
samtools view -hS -F 4 $cache/${Sample}.bt2out.sam | sed '/XS:/d' > $cache/${Sample}.single.sam
samtools view -hS -F 4 $cache/${Sample}.ecoli.bt2out.sam | sed '/XS:/d' > $cache/${Sample}.ecoli.single.sam

echo "Making bam file, sorting... "
samtools view -bS $cache/${Sample}.single.sam | samtools sort > $cache/${Sample}.sorted.bam
samtools view -bS $cache/${Sample}.ecoli.single.sam | samtools sort > $cache/${Sample}.ecoli.sorted.bam

echo "Removing duplicates... "
java -Xmx2048m -jar $PICARD/picard-2.8.0.jar MarkDuplicates \
	INPUT=$cache/${Sample}.sorted.bam \
	OUTPUT=alignments/${Sample}.uniq.bam \
	METRICS_FILE=alignments/qc/${Sample}_metrics.txt \
	REMOVE_DUPLICATES=true \
	ASSUME_SORTED=true
wait
java -Xmx2048m -jar $PICARD/picard-2.8.0.jar MarkDuplicates \
	INPUT=$cache/${Sample}.ecoli.sorted.bam \
	OUTPUT=alignments/${Sample}.ecoli.uniq.bam \
	METRICS_FILE=alignments/qc/${Sample}_ecoli_metrics.txt \
	REMOVE_DUPLICATES=true \
	ASSUME_SORTED=true
wait
	
echo "Making index file... "
samtools index alignments/${Sample}.uniq.bam
samtools index alignments/${Sample}.ecoli.uniq.bam

echo "Counting reads per chromosome..."
samtools view alignments/${Sample}.uniq.bam | cut -f 3 | sort | uniq -c | sed -e 's/^[ \t]*//' > alignments/qc/${Sample}_chromMap.txt
samtools view alignments/${Sample}.ecoli.uniq.bam | cut -f 3 | sort | uniq -c | sed -e 's/^[ \t]*//' > alignments/qc/${Sample}_ecoli_chromMap.txt

echo "Assessing fragment size distribution..."
samtools view -F 0x04 alignments/${Sample}.uniq.bam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' > alignments/qc/${Sample}_fragmentLength.txt
samtools view -F 0x04 alignments/${Sample}.ecoli.uniq.bam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' > alignments/qc/${Sample}_ecoli_fragmentLength.txt

seqDepth=$(cat alignments/qc/${Sample}.ecoli.seqDepth)
scaleFactor=`echo "1000000 / $seqDepth" | bc -l`
echo "seqDepth for " $Sample " is " $seqDepth "; scalingFactor is " $scaleFactor

# ### make bedgraph file
echo "Making bedgraph file..."
bamCoverage \
	-b alignments/${Sample}.uniq.bam \
	-o alignments/bdg/${Sample}.scaleFactor.bedgraph \
	-of bedgraph \
	-e \
	-p 4 \
	--scaleFactor $scaleFactor
	# --normalizeUsing CPM

# ### make bigwig file
echo "Making bigwig file..."
bamCoverage \
	-b alignments/${Sample}.uniq.bam \
	-o alignments/bw/${Sample}.scaleFactor.bigwig \
	-of bigwig \
	-e \
	-p 4 \
	--scaleFactor $scaleFactor
	# --normalizeUsing CPM

# ### make bedgraph file
echo "Making bedgraph file..."
bamCoverage \
	-b alignments/${Sample}.uniq.bam \
	-o alignments/bdg/${Sample}.CPM.bedgraph \
	-of bedgraph \
	-e \
	-p 4 \
	--normalizeUsing CPM

# ### make bigwig file
echo "Making bigwig file..."
bamCoverage \
	-b alignments/${Sample}.uniq.bam \
	-o alignments/bw/${Sample}.CPM.bigwig \
	-of bigwig \
	-e \
	-p 4 \
	--normalizeUsing CPM

## call peaks with macs2
macs2 callpeak -t alignments/${Sample}.uniq.bam -f BAMPE -g mm -q 0.05 -n $Sample --outdir peaks/macs/indiv/bampe
macs2 callpeak -t alignments/${Sample}.uniq.bam -f BAM -g mm -q 0.05 --shift 100 --extsize 200 -n $Sample --outdir peaks/macs/indiv/shift

### call peaks with seacr
seacr="SEACR/SEACR_1.3.sh"
controlSample=MEC_IgG_R1

bash $seacr \
	alignments/bdg/${Sample}.scaleFactor.bedgraph \
    alignments/bdg/${controlSample}.scaleFactor.bedgraph \
    non stringent \
	peaks/seacr/indiv/${Sample}_seacr_control.peaks

bash $seacr \
	alignments/bdg/${Sample}.scaleFactor.bedgraph \
	0.01 non stringent \
	peaks/seacr/indiv/${Sample}_seacr_top0.01.peaks

# Write the mapping statistics file
echo "Mapping stats..."
TOTAL=$(expr $(wc -l $cache/${Sample}_merged_R1.fastq | cut -f 1 -d' ') / 4)
FILTER=$(expr $(wc -l $cache/${Sample}_R1.trimmed.fastq | cut -f 1 -d' ') / 4)
UNMAPPED=$(samtools view -S -f 4 $cache/${Sample}.bt2out.sam | cut -f 1 | uniq | wc -l) # unmapped reads
MAPPED=$(samtools view -S -F 4 $cache/${Sample}.bt2out.sam | cut -f 1 | uniq | wc -l) # mapped reads
ONEALIGN=$(samtools view -S $cache/${Sample}.single.sam | grep -vc XS: ) #  reads mapping to one postion
UNIQ=$(samtools view alignments/${Sample}.uniq.bam | wc -l) # unique reads

printf "TOTAL\tFILTER\tUNMAPPED\tMAPPED\tONEALIGN\tUNIQ\n" >  alignments/qc/${Sample}_map_stats.txt
printf "%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\n" $TOTAL $FILTER $UNMAPPED $MAPPED $ONEALIGN $UNIQ >> alignments/qc/${Sample}_map_stats.txt

echo "Mapping ecoli stats..."
TOTAL=$(expr $(wc -l $cache/${Sample}_merged_R1.fastq | cut -f 1 -d' ') / 4)
FILTER=$(expr $(wc -l $cache/${Sample}_R1.trimmed.fastq | cut -f 1 -d' ') / 4)
UNMAPPED=$(samtools view -S -f 4 $cache/${Sample}.ecoli.bt2out.sam | cut -f 1 | uniq | wc -l) # unmapped reads
MAPPED=$(samtools view -S -F 4 $cache/${Sample}.ecoli.bt2out.sam | cut -f 1 | uniq | wc -l) # mapped reads
ONEALIGN=$(samtools view -S $cache/${Sample}.ecoli.single.sam | grep -vc XS: ) #  reads mapping to one postion
UNIQ=$(samtools view alignments/${Sample}.ecoli.uniq.bam | wc -l) # unique reads

printf "TOTAL\tFILTER\tUNMAPPED\tMAPPED\tONEALIGN\tUNIQ\n" >  alignments/qc/${Sample}_ecoli_map_stats.txt
printf "%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\n" $TOTAL $FILTER $UNMAPPED $MAPPED $ONEALIGN $UNIQ >> alignments/qc/${Sample}_ecoli_map_stats.txt


echo "Doing MultiQC... "
multiqc .

echo "Cleaning up... "
rm $cache/${Sample}_1.fastq.gz
rm $cache/${Sample}_2.fastq.gz
rm $cache/${Sample}_1.trimmed.fastq.gz
rm $cache/${Sample}_2.trimmed.fastq.gz
rm $cache/${Sample}_1.unpaired.fastq.gz
rm $cache/${Sample}_2.unpaired.fastq.gz
rm $cache/${Sample}.bt2out.sam
rm $cache/${Sample}.single.sam
rm $cache/${Sample}.sorted.bam

mv *.out logs
mv *.err logs

echo "Done!"
