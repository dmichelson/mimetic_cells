#!/bin/bash

#SBATCH -c 4
#SBATCH -N 1
#SBATCH -t 0-36:00
#SBATCH -p medium
#SBATCH --mem-per-cpu=16G
#SBATCH -o logs/dam41_%j.out
#SBATCH -e logs/dam41_%j.err

module load gcc/6.2.0
#module load python/2.7.12
#module load java/jdk-1.8u112
#module load git/2.9.5
#module load bwa/0.7.8
#module load cuda/9.0
#module load sratoolkit/2.8.1
#module load fastqc/0.11.5
module load bowtie2/2.2.9
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

cd /n/groups/cbdm-db/dam41/mec_scatac

source activate snaptools_env

prefix=$1

snaptools dex-fastq \
	--input-fastq=data/200305_A00794_0197_BHHM5WDRXX/fastq/mTEC_scATAC_1/$prefix/${prefix}_L001_R1_001.fastq.gz \
	--output-fastq=data/dex/${prefix}_L001_R1.dex.fastq.gz \
	--index-fastq-list data/200305_A00794_0197_BHHM5WDRXX/fastq/mTEC_scATAC_1/$prefix/${prefix}_L001_R2_001.fastq.gz

snaptools dex-fastq \
	--input-fastq=data/200305_A00794_0197_BHHM5WDRXX/fastq/mTEC_scATAC_1/$prefix/${prefix}_L001_R3_001.fastq.gz \
	--output-fastq=data/dex/${prefix}_L001_R3.dex.fastq.gz \
	--index-fastq-list data/200305_A00794_0197_BHHM5WDRXX/fastq/mTEC_scATAC_1/$prefix/${prefix}_L001_R2_001.fastq.gz

snaptools dex-fastq \
	--input-fastq=data/200305_A00794_0197_BHHM5WDRXX/fastq/mTEC_scATAC_1/$prefix/${prefix}_L002_R1_001.fastq.gz \
	--output-fastq=data/dex/${prefix}_L002_R1.dex.fastq.gz \
	--index-fastq-list data/200305_A00794_0197_BHHM5WDRXX/fastq/mTEC_scATAC_1/$prefix/${prefix}_L002_R2_001.fastq.gz

snaptools dex-fastq \
	--input-fastq=data/200305_A00794_0197_BHHM5WDRXX/fastq/mTEC_scATAC_1/$prefix/${prefix}_L002_R3_001.fastq.gz \
	--output-fastq=data/dex/${prefix}_L002_R3.dex.fastq.gz \
	--index-fastq-list data/200305_A00794_0197_BHHM5WDRXX/fastq/mTEC_scATAC_1/$prefix/${prefix}_L002_R2_001.fastq.gz

cat data/dex/${prefix}_L001_R1.dex.fastq.gz data/dex/${prefix}_L002_R1.dex.fastq.gz > data/dex/${prefix}_R1.dex.fastq.gz	
cat data/dex/${prefix}_L001_R3.dex.fastq.gz data/dex/${prefix}_L002_R3.dex.fastq.gz > data/dex/${prefix}_R3.dex.fastq.gz	

genomeDir="/n/groups/shared_databases/bowtie2_indexes"
bowtie2Dir="/n/app/bowtie2/2.2.9/bin/"
genome="mm10"

snaptools align-paired-end \
	--input-reference=$genomeDir/$genome \
	--input-fastq1=data/dex/${prefix}_R1.dex.fastq.gz \
	--input-fastq2=data/dex/${prefix}_R3.dex.fastq.gz \
	--output-bam=alignments/${prefix}.bam \
	--aligner=bowtie2 \
	--path-to-aligner=$bowtie2Dir \
	--read-fastq-command=zcat \
	--min-cov=0 \
	--num-threads=4 \
	--if-sort=True \
	--tmp-folder=./temp \
	--overwrite=TRUE
	
fetchChromSizes $genome > $genome.gs

snaptools snap-pre  \
	--input-file=alignments/${prefix}.bam  \
	--output-snap=alignments/${prefix}.snap  \
	--genome-name=$genome  \
	--genome-size=reference/$genome.chrom.sizes  \
	--min-mapq=30  \
	--min-flen=0  \
	--max-flen=1000  \
	--keep-chrm=TRUE  \
	--keep-single=FALSE  \
	--keep-secondary=FALSE  \
	--overwrite=True  \
	--min-cov=100  \
	--verbose=True
	
snaptools snap-add-bmat	\
	--snap-file=alignments/${prefix}.snap	\
	--bin-size-list 1000 5000 10000	\
	--verbose=True

snaptools snap-add-pmat \
	--snap-file=alignments/${prefix}.snap \
	--peak-file=peaks/peaks.combined.bed