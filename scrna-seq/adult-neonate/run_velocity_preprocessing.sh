#!/bin/bash

#SBATCH -c 8
#SBATCH -N 1
#SBATCH -t 0-2:00
#SBATCH -p short
#SBATCH --mem-per-cpu=4G
#SBATCH -o logs/dam41_%j.out
#SBATCH -e logs/dam41_%j.err

cd /n/groups/cbdm-db/your/working/directory

## create conda env and install dependencies
conda create -n rna-velocity python numpy scipy cython numba matplotlib scikit-learn h5py click
pip install git+https://github.com/pachterlab/kb_python@devel

source activate rna-velocity
# conda deactivate rna-velocity

## make kallisto-bustools reference
mkdir kb
mkdir kb-count-out
mkdir scvelo
mkdir seurat

cd kb
wget ftp://ftp.ensembl.org/pub/release-101/gtf/mus_musculus/Mus_musculus.GRCm38.101.gtf.gz
wholeGenomeFasta="/n/groups/shared_databases/igenome/Mus_musculus/NCBI/GRCm38/Sequence/WholeGenomeFasta/genome.fa"
kb ref -i index.idx -g t2g.txt \
    -f1 cdna.fa -f2 intron.fa \
    -c1 cdna_t2c.txt -c2 intron_t2c.txt \
    --workflow lamanno -n 4 \
    --tmp /n/scratch3/users/YOUR/SCRATCH/FOLDER \
    $wholeGenomeFasta Mus_musculus.GRCm38.101.gtf

## make loom file with kallisto-bustools
cd ..
fastqDir="/n/groups/cbdm_lab/YOUR/FASTQ/FILES"
kb count -i kb/index.idx_cdna,kb/index.idx_intron.0,kb/index.idx_intron.1,kb/index.idx_intron.2 \
    -g kb/t2g.txt \
    -x 10xv3 --workflow lamanno --h5ad -m 30G -t 8 \
    --tmp /n/scratch3/users/YOUR/SCRATCH/FOLDER \
    -c1 kb/cdna_t2c.txt -c2 kb/intron_t2c.txt \
    -o kb-count-out \
    --overwrite \
    $fastqDir/YOUR_FASTQ_FILE_R1_001.fastq.gz $fastqDir/YOUR_FASTQ_FILE_R2_001.fastq.gz

